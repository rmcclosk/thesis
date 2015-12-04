#!/usr/bin/env python3

import copy
import hashlib
import itertools
import logging
import os
import pprint
import random
import subprocess
import sqlite3
import sys
import tempfile
import time
import yaml

class Tempfile:
    """A temporary file created with mkstemp()"""
    def __init__(self, mode="r"):
        self.fd, self.name = tempfile.mkstemp()
        self.handle = os.fdopen(self.fd, mode)

    def write(self, text):
        self.handle.write(text)

    def close(self):
        self.handle.close()

    def delete(self):
        os.unlink(self.name)


class ParameterSet:
    """A set of experimental parameters to be combinatorially expanded"""
    def __init__(self, spec={}, exclusions=[]):
        self.parameters = spec

        # parse things like primitives and "range(10)"
        for param, value in self.parameters.items():
            if isinstance(value, str):
                try:
                    value = eval(value)
                except (NameError, SyntaxError):
                    pass
            if isinstance(value, str) or not hasattr(value, '__iter__'):
                value = [value]
            if not isinstance(value, str) and hasattr(value, '__iter__'):
                value = list(value)
            value = list(map(str, value))
            self.parameters[param] = value

        self.expand = dict.fromkeys(self.parameters.keys(), True)
        self.exclusions = [ParameterSet(excl) for excl in exclusions]

    def __len__(self):
        return len(self.parameters)

    def __iter__(self):
        expand_params = {k: v for k, v in self.parameters.items() if self.expand[k]}
        no_expand_params = {k: v for k, v in self.parameters.items() if not self.expand[k]}
        for values in itertools.product(*expand_params.values()):
            p = dict(zip(expand_params.keys(), ([v] for v in values)))
            p.update(no_expand_params)
            if not self.is_excluded(p):
                yield p

    def __getitem__(self, key):
        return self.parameters[key]

    def __setitem__(self, key, value):
        self.parameters[key] = value

    def __str__(self):
        expand_params = {k: v for k, v in self.parameters.items() if self.expand[k]}
        return pprint.pformat(expand_params, width=120, compact=True)

    @property
    def names(self):
        return self.parameters.keys()

    def __iadd__(self, other):
        my_names = set(self.names)
        other_names = set(other.names)

        # exclusions are carried forward
        self.exclusions.extend(other.exclusions)

        # all pararameters from the other set, which are not in this set, get
        # copied but not expanded
        for name in other_names - my_names:
            self[name] = copy.deepcopy(other[name])
            self.expand[name] = False

        # for shared names
        for name in my_names & other_names:
            if not set(self[name]).issubset(set(other[name])):
                raise ValueError("Cannot combine parameter {}".format(name))

        return self

    def __and__(self, other):
        new_set = ParameterSet({})
        for key in self.names:
            new_set[key] = list(set(self[key]) & set(other[key]))
            new_set.expand[key] = self.expand[key]
        return new_set

    def __eq__(self, other):
        if len(self) != len(other):
            return False
        for key in self.names:
            if not key in other.names:
                return False
            if len(set(self[key]) ^ set(other[key])) > 0:
                return False
        return True

    def __contains__(self, params):
        return (params & self == params) and not self.is_excluded(params.parameters)

    def is_excluded(self, params):
        for excl in self.exclusions:
            for e in excl:
                if ParameterSet(e) in ParameterSet(params):
                    return True
        return False


class Step:
    """One step of a computational experiment."""
    def __init__(self, name, spec, globals, cur):
        """Initialize a step from a YAML spec"""
        self.name = name
        self.dir = os.path.join(globals["Name"], self.name)
        try:
            os.makedirs(self.dir)
        except FileExistsError:
            pass
        self.cur = cur

        self.extension = spec["Extension"]
        self.interpreter = spec["Interpreter"]
        self.startup = spec.get("Startup", "")
        self.rule = spec["Rule"]
        self.sleep = spec.get("Sleep", globals.get("Sleep", 60))
        self.nproc = spec.get("Processes", globals.get("Processes", 1))

        self.parameters = ParameterSet(spec.get("Parameters", {}),
                                       spec.get("Exclusions", []))

    def set_dependencies(self, depends):
        self.depends = depends

    def find_file(self, parameters):
        """Find the name for a file with a particular set of parameters."""
        query = "SELECT path FROM data WHERE step = ? AND parameter = ? AND value = ?"
        query = " INTERSECT ".join([query] * len(parameters))
        query_params = []
        for k, v in parameters.items():
            query_params.extend([self.name, k, str(v)])
        self.cur.execute(query, tuple(query_params))

        result = self.cur.fetchone()
        if result is not None:
            return result[0]

        basename = hashlib.md5(bytes(str(parameters), "UTF-8")).hexdigest()
        filename = "{}.{}".format(basename, self.extension)
        return os.path.join(self.dir, filename)

    def expand_parameters(self, parameters):
        """Add file names and extras to parameters"""
        for k, v in parameters.items():
            parameters[k] = " ".join(v)
        parameters["yaml"] = yaml.dump(parameters, width=1000).rstrip()
        parameters["seed"] = random.randrange(2**31)

    def store_parameters(self, path, parameters):
        """Add a path to the database associated with some parameters."""
        for k, v in parameters.items():
            self.cur.execute("INSERT OR REPLACE INTO data (step, path, parameter, value) VALUES (?, ?, ?, ?)",
                             (self.name, path, k, str(v)))

    def store_checksum(self, path):
        """Store the MD5 checksum of a path in the database"""
        if not os.path.exists(path):
            logging.warning("Path {} was not created".format(path))
            return

        with open(path, "rb") as f:
            checksum = hashlib.md5(f.read()).hexdigest()
            self.cur.execute("INSERT OR REPLACE INTO md5 (path, checksum) VALUES (?, ?)",
                             (path, checksum))

    def store_dependencies(self, path, depends):
        """Store checksums for dependencies used to make a file."""
        for dep_path in itertools.chain(*depends.values()):
            self.cur.execute("SELECT checksum FROM md5 WHERE path = ?", (dep_path,))
            self.cur.execute("""
                INSERT OR REPLACE INTO dependencies (path, dep_path, dep_checksum)
                VALUES (?, ?, ?)
            """, (path, dep_path, self.cur.fetchone()[0]))

    def get_prerequisites(self, parameters):
        """Find all the prerequisite files to build a target."""
        prereqs = {}
        for dep_step in self.depends:
            prereqs[dep_step.name] = []
            for p in dep_step.parameters & parameters:
                matching_file = dep_step.find_file(p)
                if not os.path.exists(matching_file):
                    raise RuntimeError("Prerequisites from step {} for step {} with parameters {} are missing"
                                       .format(self.name, dep_step.name, p))
                prereqs[dep_step.name].append(matching_file)
        return prereqs

    def needs_update(self, path, depends):
        """Check if a target needs to be remade."""
        self.cur.execute("SELECT checksum FROM md5 WHERE path = ?", (path,))
        if self.cur.fetchone() is None:
            logging.info("Target {} does not exist".format(path))
            return True

        for dep_path in itertools.chain(*depends.values()):
            self.cur.execute("SELECT checksum == dep_checksum FROM md5 JOIN dependencies ON md5.path == dependencies.dep_path WHERE dependencies.path = ? AND dep_path = ?",
                             (path, dep_path))
            if not self.cur.fetchone()[0]:
                logging.info("Dependency {} for target {} has changed".format(dep_path, path))
                return True

        logging.info("File {} is up to date".format(path))
        return False

    def run(self):
        """Run a step."""
        targets = []

        scripts = []
        script_cur = 0

        # add a copy of the rule for each set of parameters, divided amongst
        # scripts for each process
        for p in self.parameters:
            target = self.find_file(p)
            prereqs = self.get_prerequisites(p)

            if self.needs_update(target, prereqs):

                if len(scripts) == script_cur:
                    scripts.append(Tempfile("w"))
                    scripts[script_cur].write("{}\n".format(self.startup))
                    targets.append([])

                self.store_parameters(target, p)
                self.store_dependencies(target, prereqs)
                targets[script_cur].append(target)

                self.expand_parameters(p)
                p[self.name] = target
                p.update({k: " ".join(v) for k, v in prereqs.items()})

                scripts[script_cur].write("{}\n".format(self.rule.format(**p)))
                logging.debug(self.rule.format(**p))

                script_cur += 1
                if (self.nproc > 0):
                    script_cur %= self.nproc

        # make a driver for each script
        drivers = []
        for script in scripts: 
            script.close()

            driver = Tempfile("w")
            driver.write("#!{}\n".format(os.environ["SHELL"]))
            driver.write("{} < {}\n".format(self.interpreter, script.name))
            driver.close()
            drivers.append(driver)

        # submit each driver
        procs = []
        for driver in drivers:
            procs.append(subprocess.Popen([os.environ["SHELL"], driver.name]))

        # wait for them to finish
        done = [False] * len(procs)
        target_cur = [0] * len(procs)

        while not all(done):
            done = [proc.poll() is not None for proc in procs]

            for i, files in enumerate(targets):
                if done[i] and proc[i].returncode == 0:
                    while target_cur[i] < len(files):
                        self.store_checksum(files[target_cur[i]])
                        target_cur[i] += 1
                    logging.info("Process {} is finished".format(files[target_cur[i]]))
                elif target_cur[i] < len(files) - 1 and os.path.exists(files[target_cur[i]+1]):
                    logging.info("File {} was created".format(files[target_cur[i]]))
                    self.store_checksum(files[target_cur[i]])
                    target_cur[i] += 1

            if not all(done):
                time.sleep(self.sleep)

        for f in scripts + drivers: f.delete()

        for proc in procs:
            if proc.returncode != 0:
                return False
        return True


class Experiment:
    """A computational experiment."""
    def __init__(self, spec):
        """Initialize an experiment from a YAML spec"""
        self.name = spec["Name"]
        self.nproc = spec.get("Processes", 1)
        self.dir = self.name
        try:
            os.makedirs(self.dir)
        except FileExistsError:
            pass

        self.con = sqlite3.connect("{}.sqlite".format(self.name))
        self.con.isolation_level = None
        self.cur = self.con.cursor()
        self.create_database()

        self.steps = {}
        self.step_graph = {}
        for step_name in spec.get("Steps", []):
            step_spec = spec["Steps"][step_name]
            self.steps[step_name] = Step(step_name, step_spec, spec, self.cur)

            depends = step_spec.get("Depends", [])
            if isinstance(depends, str):
                depends = depends.split(" ")
            self.step_graph[step_name] = depends

        # now that we've read them all, we can set all the steps' dependencies
        for step_name in self.step_graph:
            depends = [self[dep_step] for dep_step in self.step_graph[step_name]]
            self[step_name].set_dependencies(depends)

        self.propagate_parameters()
        self.sanitize_database()
        self.sanitize_filesystem()

    def __iter__(self):
        dag = copy.deepcopy(self.step_graph)
        for i in range(len(dag)):
            next_step = min(dag.items(), key = lambda x: len(x[1]))[0]
            for step in dag:
                try:
                    dag[step].remove(next_step)
                except ValueError:
                    pass
            del dag[next_step]
            yield self.steps[next_step]

    def __getitem__(self, key):
        return self.steps[key]

    def propagate_parameters(self):
        """Propagate steps' prerequisites' parameters."""
        for step in self:
            for dep_step in self.step_graph[step.name]:
                step.parameters += self[dep_step].parameters

    def create_database(self):
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS md5 (
                path TEXT PRIMARY KEY,
                checksum TEXT
            )
        """)
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS data (
                step TEXT,
                path TEXT,
                parameter TEXT, 
                value TEXT,
                PRIMARY KEY (step, path, parameter))
        """)
        self.cur.execute("""
            CREATE TABLE IF NOT EXISTS dependencies (
                path TEXT,
                dep_path TEXT,
                dep_checksum TEXT,
                PRIMARY KEY (path, dep_path))
        """)

    def purge(self, path):
        """Delete a path from the database."""
        self.cur.execute("DELETE from md5 WHERE path = ?", (path,))
        self.cur.execute("DELETE from data WHERE path = ?", (path,))
        self.cur.execute("DELETE from dependencies WHERE path = ?", (path,))

        # delete anybody that depended on this file too
        self.cur.execute("SELECT path FROM dependencies WHERE dep_path = ?", (path,))
        other_paths = [row[0] for row in self.cur.fetchall()]
        for p in other_paths:
            logging.info("Recursively deleting {}".format(p))
            self.purge(p)

    def sanitize_database(self):
        """Clean up the database."""

        # make sure all the tables are in sync
        self.cur.execute("SELECT path FROM md5")
        md5_paths = set([row[0] for row in self.cur.fetchall()])
    
        self.cur.execute("SELECT DISTINCT path FROM data")
        data_paths = set([row[0] for row in self.cur.fetchall()])
    
        self.cur.execute("SELECT DISTINCT path FROM dependencies")
        dep_paths = set([row[0] for row in self.cur.fetchall()])
    
        for path in md5_paths | data_paths | dep_paths:
            if not os.path.exists(path):
                logging.info("Deleting file {} from the database (not on filesystem)".format(path))
                self.purge(path)
            else:
                with open(path, "rb") as f:
                    checksum = hashlib.md5(f.read()).hexdigest()
                self.cur.execute("SELECT checksum FROM md5 WHERE path = ?", (path,))
                row = self.cur.fetchone()
                if row is None or checksum != row[0]:
                    logging.info("Deleting file {} from the database (checksum doesn't match, or no checksum in database)".format(path))
                    self.purge(path)
    
        for path in (md5_paths ^ data_paths) | (dep_paths - md5_paths - data_paths):
            logging.info("Deleting file {} from the database (tables out of sync)".format(path))
            self.purge(path)

        for step in self:
            # get all paths in the database from this step
            self.cur.execute("SELECT DISTINCT path FROM data WHERE step = ?", (step.name,))
            paths = set([row[0] for row in self.cur.fetchall()])

            # for each path, see if there is a matching set of parameters
            # if not, delete the path
            for path in paths:
                self.cur.execute("SELECT parameter, value FROM data WHERE path = ?", (path,))
                params = ParameterSet({row[0]: eval(row[1]) for row in self.cur.fetchall()})
                if not params in step.parameters:
                    logging.info("Deleting file {} from the database (unused parameters)".format(path))
                    self.purge(path)

    def sanitize_filesystem(self):
        """Clean the file system."""
        for dirpath, dirnames, filenames in os.walk(self.dir):
            for basename in filenames:
                path = os.path.join(dirpath, basename)
                self.cur.execute("SELECT * FROM md5 WHERE path = ?", (path,))
                if self.cur.fetchone() is None:
                    logging.info("Deleting file {} from filesystem (not in database)".format(path))
                    os.unlink(path)

    def run(self):
        """Run an experiment."""
        for step in self:
            logging.info("Starting step {}".format(step.name))
            if not step.run():
                logging.error("Step {} failed".format(step.name))
                break


with open("test.yaml") as f:
    logging.basicConfig(level=logging.DEBUG)
    e = Experiment(yaml.load(f))
    e.run()
