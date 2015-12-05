#!/usr/bin/env python3

import copy
import itertools
import os
import pprint
import yaml

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
                yield {k: v.pop() for k, v in d.items()}

    def filter(self, other):
        """TODO"""
        pass

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

        self.parameters = ParameterSet(spec.get("Parameters", {}))
        self.exclusions = [ParameterSet(excl) for excl in spec.get("Exclusions", [])]
        self.depends = [] # must be filled in later

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


with open("test.yaml") as f:
    expt = yaml.load(f)
    spec = expt["Steps"]["network"]
    s = Step("network", spec, expt, 0)
    for p in s.parameters.combinations(s.exclusions):
        print(p)
