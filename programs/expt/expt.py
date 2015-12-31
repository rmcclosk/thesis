#!/usr/bin/env python3

import argparse
import copy
import datetime
import hashlib
import itertools
import logging
import os
import pprint
import random
import re
import socket
import sqlite3
import subprocess
import sys
import tempfile
import time
import yaml

class Tempfile:
    """A temporary file created with mkstemp()"""
    def __init__(self, mode="r", tmpdir=None):
        self.fd, self.name = tempfile.mkstemp(dir=tmpdir)
        self.handle = os.fdopen(self.fd, mode)

    def write(self, text):
        self.handle.write(text)

    def close(self):
        self.handle.close()

    def delete(self):
        os.unlink(self.name)


class ParameterSet (dict):
    """A set of experimental parameters to be combinatorially expanded"""
    def __init__(self, spec={}):
        """Create a ParameterSet from a dictionary specification"""
        for param, value in spec.items():
            if isinstance(value, str):
                try:
                    value = eval(value)
                except (NameError, SyntaxError):
                    pass
            if isinstance(value, str) or not hasattr(value, '__iter__'):
                value = [value]
            self[param] = set(value)

    def intersect(self, other):
        """Get the parameters and values which are in common with myself and another set"""
        new_set = ParameterSet()
        for name in self:
            if name in other and (self[name] & other[name]):
                new_set[name] = self[name] & other[name]
        return new_set

    def combinations(self, exclusions = []):
        """Iterate over all possible combinations of parameters"""
        for values in itertools.product(*self.values()):
            d = ParameterSet(dict(zip(self.keys(), values)))

            # check if this one is excluded
            keep = True
            for excl in exclusions:
                if len(d.intersect(excl)) == len(excl):
                    keep = False
                    break

            if keep: 
                yield ParameterSet({k: v.pop() for k, v in d.items()})

    def filter(self, other):
        """Find all the parameters in one set which match another"""
        new_set = copy.deepcopy(self)
        for key in other:
            if key in self:
                new_set[key] &= other[key]
        return new_set

    def as_dict(self):
        d = {}
        for key in self:
            if len(self[key]) == 1:
                d[key] = list(self[key])[0]
            else:
                d[key] = " ".join(sorted(map(str, self[key])))
        return d


class Task:
    """A task running either in a shell, or on a cluster"""
    def __init__(self, driver, qsub=False):
        self.qsub = qsub
        if qsub:
            proc = subprocess.Popen(["qsub", "-V", driver.name], stdout=subprocess.PIPE)
            self.id = str(proc.communicate()[0], "UTF-8").split(".")[0]
        else:
            self.proc = subprocess.Popen([os.environ["SHELL"], driver.name])

    def finished(self):
        if self.qsub:
            proc = subprocess.Popen(["qstat", "-u", os.environ["USER"]], stdout=subprocess.PIPE)
            out = str(proc.communicate()[0], "UTF-8")
            for line in out.split("\n"):
                line = line.split()
                if len(line) > 9 and line[9] != "C" and self.id in line[0]:
                    return False
            return True
        else:
            return self.proc.poll() is not None

    def returncode(self):
        if self.qsub:
            return 0
        return self.proc.returncode


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
        self.nthread = spec.get("Threads", globals.get("Threads", 1))
        self.memory = spec.get("Memory", globals.get("Memory", "512m"))
        walltime = spec.get("Walltime", globals.get("Walltime", "00:05:00"))
        try:
            h, m, s = map(int, walltime.split(":"))
        except ValueError:
            logging.error("Improperly formatted walltime {}".format(walltime))
            raise
        self.walltime = datetime.timedelta(hours=h, minutes=m, seconds=s)

        self.parameters = ParameterSet(spec.get("Parameters", {}))
        self.exclusions = [ParameterSet(excl) for excl in spec.get("Exclusions", [])]

        self.found_params = []
        self.found_files = []

    def __contains__(self, parameters):
        """Check if a ParameterSet will be used by this step"""
        if set(parameters.keys()) != set(self.parameters.keys()):
            return False
        for key in parameters:
            if not parameters[key].issubset(self.parameters[key]):
                return False
        for excl in self.exclusions:
            if len(parameters.intersect(excl)) == len(excl):
                return False
        return True

    def get_walltime(self, ntasks):
        """Get a string specifying the amount of walltime for n tasks"""
        wt = self.walltime * ntasks
        h = int(wt.seconds / 3600)
        m = int(wt.seconds % 3600 / 60)
        s = wt.seconds % 60
        return "{:02d}:{:02d}:{:02d}".format(h, m, s)

    def find_file(self, parameters):
        """Find the name for a file with a particular set of parameters."""
        pdict = parameters.as_dict()
        try:
            idx = self.found_params.index(pdict)
            return self.found_files[idx]
        except ValueError:
            pass

        query = "SELECT path FROM data WHERE step = ? AND parameter = ? AND value = ?"
        query = " INTERSECT ".join([query] * len(parameters))
        query_params = []
        for k, v in parameters.as_dict().items():
            query_params.extend([self.name, k, str(v)])
        self.cur.execute(query, tuple(query_params))

        result = self.cur.fetchone()
        if result is not None:
            self.found_params.append(pdict)
            self.found_files.append(result[0])
            logging.debug("Matching file found in database for parameters {}".format(parameters.as_dict()))
            return result[0]

        logging.debug("No matching file found in database for parameters {}".format(parameters.as_dict()))
        basename = hashlib.md5(bytes(str(parameters), "UTF-8")).hexdigest()
        filename = "{}.{}".format(basename, self.extension)
        self.found_params.append(pdict)
        self.found_files.append(os.path.join(self.dir, filename))
        return os.path.join(self.dir, filename)

    def store_parameters(self, path, parameters):
        """Add a path to the database associated with some parameters."""
        for k, v in parameters.as_dict().items():
            self.cur.execute("INSERT OR REPLACE INTO data (step, path, parameter, value) VALUES (?, ?, ?, ?)",
                             (self.name, path, k, str(v)))

    def store_checksum(self, path):
        """Store the MD5 checksum of a path in the database"""
        if not os.path.exists(path):
            logging.warning("Path {} was not created".format(path))
            return

        with open(path, "rb") as f:
            checksum = hashlib.md5(f.read()).hexdigest()
            logging.debug("Storing checksum for file {}".format(path))
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

    def get_prerequisites(self, parameters, dep_steps):
        """Find all the prerequisite files to build a target."""
        prereqs = {}
        for dep_step in dep_steps:
            prereqs[dep_step.name] = []
            dep_params = dep_step.parameters.filter(parameters)
            dep_excl = self.exclusions + dep_step.exclusions
            for p in dep_params.combinations(dep_excl):
                matching_file = dep_step.find_file(p)
                if not os.path.exists(matching_file):
                    logging.warning("Prerequisites from step {} for step {} with parameters {} are missing"
                                    .format(self.name, dep_step.name, p.as_dict()))
                else:
                    prereqs[dep_step.name].append(matching_file)
            prereqs[dep_step.name] = sorted(prereqs[dep_step.name])
            if len(prereqs[dep_step.name]) == 0:
                logging.info("Parameters {} for step {} do not match any files from step {}"
                             .format(parameters.as_dict(), self.name, dep_step.name))
        return prereqs

    def needs_update(self, path, depends):
        """Check if a target needs to be remade."""
        if any(len(x) == 0 for x in depends.values()):
            logging.info("No prerequisites exist for target {}".format(path))
            return False

        self.cur.execute("SELECT checksum FROM md5 WHERE path = ?", (path,))
        if self.cur.fetchone() is None:
            logging.info("Target {} does not exist".format(path))
            return True

        for dep_path in itertools.chain(*depends.values()):
            self.cur.execute("""
                SELECT checksum == dep_checksum FROM md5
                JOIN dependencies ON md5.path == dependencies.dep_path 
                WHERE dependencies.path = ? AND dep_path = ?""",
                (path, dep_path))
            row = self.cur.fetchone()
            if row is None or not row[0]:
                logging.info("Dependency {} for target {} has changed".format(dep_path, path))
                return True

        logging.info("File {} is up to date".format(path))
        return False

    def run(self, dep_steps=[], qsub=False, tmpdir=None):
        """Run a step."""

        targets = []
        scripts = []
        script_cur = 0

        # add a copy of the rule for each set of parameters, divided amongst
        # scripts for each process
        for p in self.parameters.combinations(self.exclusions):
            target = self.find_file(p)
            prereqs = self.get_prerequisites(p, dep_steps)

            if self.needs_update(target, prereqs):

                if len(scripts) == script_cur:
                    scripts.append(Tempfile("w", tmpdir))
                    scripts[script_cur].write("{}\n".format(self.startup))
                    targets.append([])

                self.store_parameters(target, p)
                self.store_dependencies(target, prereqs)
                targets[script_cur].append(target)

                p = p.as_dict()
                p["yaml"] = yaml.dump(p, width=1000).rstrip()
                p["seed"] = random.randrange(2**31)
                p["$#"] = sum(len(x) for x in prereqs.values())
                p[self.name] = target
                p.update({k: " ".join(sorted(v)) for k, v in prereqs.items()})

                scripts[script_cur].write("{}\n".format(self.rule.format(**p)))
                logging.debug(self.rule.format(**p))

                script_cur += 1
                if (self.nproc > 0):
                    script_cur %= self.nproc

        # make a driver for each script
        drivers = []
        for i, script in enumerate(scripts):
            script.close()

            driver = Tempfile("w", tmpdir)
            driver.write("#!{}\n".format(os.environ["SHELL"]))
            if qsub:
                driver.write("#PBS -S {}\n".format(os.environ["SHELL"]))
                driver.write("#PBS -l nodes=1:ppn={}\n".format(self.nthread))
                walltime = self.get_walltime(len(targets[i]))
                driver.write("#PBS -l walltime={}\n".format(walltime))
                driver.write("#PBS -l pmem={}\n".format(self.memory))
                logging.debug("Submitted job {} with ppn={}, walltime={}, pmem={}"
                              .format(i, self.nthread, walltime, self.memory))
                driver.write("cd $PBS_O_WORKDIR\n")
            driver.write("{} < {}\n".format(self.interpreter, script.name))
            driver.close()
            drivers.append(driver)

        # submit each driver
        tasks = []
        for driver in drivers:
            tasks.append(Task(driver, qsub))

        # wait for them to finish
        done = [False] * len(tasks)
        target_cur = [0] * len(tasks)

        while not all(done):
            done = [task.finished() for task in tasks]

            for i, files in enumerate(targets):
                if done[i] and tasks[i].returncode() == 0:
                    while target_cur[i] < len(files):
                        self.store_checksum(files[target_cur[i]])
                        target_cur[i] += 1
                    logging.info("Process {} is finished".format(i))
                else:
                    while target_cur[i] < len(files) - 1 and os.path.exists(files[target_cur[i]+1]):
                        logging.info("File {} was created".format(files[target_cur[i]]))
                        self.store_checksum(files[target_cur[i]])
                        target_cur[i] += 1

            if not all(done):
                time.sleep(self.sleep)

        for f in scripts + drivers: f.delete()

        for task in tasks:
            if task.returncode() != 0:
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
        for step_name in spec.get("Steps", []):
            step_spec = spec["Steps"][step_name]
            self.steps[step_name] = Step(step_name, step_spec, spec, self.cur)

        self.step_graph = {}
        for step_name in spec.get("Steps", []):
            step_spec = spec["Steps"][step_name]
            depends = step_spec.get("Depends", [])
            if isinstance(depends, str):
                depends = depends.split(" ")
            self.step_graph[step_name] = depends

        self.sanitize_database()
        self.sanitize_filesystem()

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
        self.cur.execute("""
            CREATE INDEX IF NOT EXISTS data_idx
            ON data (step, path, parameter, value)
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

        for step in self.steps.values():
            # get all paths in the database from this step
            self.cur.execute("SELECT DISTINCT path FROM data WHERE step = ?", (step.name,))
            paths = set([row[0] for row in self.cur.fetchall()])

            # for each path, see if there is a matching set of parameters
            # if not, delete the path
            for path in paths:
                self.cur.execute("SELECT parameter, value FROM data WHERE path = ?", (path,))
                params = ParameterSet({row[0]: row[1] for row in self.cur.fetchall()})
                for p in params:
                    try:
                        params[p] = eval(str(params[p]))
                    except (NameError, SyntaxError):
                        pass
                if not params in step:
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

    def iter_steps(self):
        dag = copy.deepcopy(self.step_graph)
        nstep = len(dag)
        for i in range(nstep):
            next_step = min(dag.items(), key = lambda x: len(x[1]))[0]
            for step in dag:
                try:
                    dag[step].remove(next_step)
                except ValueError:
                    pass
            del dag[next_step]
            yield self.steps[next_step]

    def run(self, qsub=False, tmpdir=None, which_step=None):
        """Run an experiment."""
        for step in self.iter_steps():
            if which_step is not None and which_step != step.name: continue
            logging.info("Starting step {}".format(step.name))
            dep_steps = [self.steps[s] for s in self.step_graph[step.name]]
            if not step.run(dep_steps, qsub, tmpdir):
                logging.error("Step {} failed".format(step.name))
                break
            logging.info("Finished step {}".format(step.name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a computational experiment.")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("-q", "--qsub", action="store_true", default=False)
    parser.add_argument("-t", "--tmpdir", default=None)
    parser.add_argument("-s", "--step")
    parser.add_argument("yaml_file", type=argparse.FileType("r"))
    args = parser.parse_args()

    loglevels = [logging.WARNING, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    spec = yaml.load(args.yaml_file)
    
    try:
        hostname_re = re.compile(spec["Hostname"])
    except KeyError:
        hostname_re = re.compile(".*")
    if hostname_re.match(socket.gethostname()) is None:
        logging.error("Wrong hostname: expected to be running on {}, but we are on {}".format(hostname_re.pattern, socket.gethostname()))
        sys.exit()

    expt = Experiment(spec)
    expt.run(args.qsub, args.tmpdir, args.step)
