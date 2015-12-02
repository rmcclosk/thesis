#!/usr/bin/env python3

import argparse
import subprocess
import datetime
import hashlib
import itertools
import logging
import os
import random
import sqlite3
import sys
import tempfile
import time
import yaml
from pprint import pprint

class Tempfile:
    """A temporary file created with mkstemp()"""
    def __init__(self, tmpdir, mode):
        self.fd, self.name = tempfile.mkstemp(dir=tmpdir)
        self.handle = os.fdopen(self.fd, mode)

    def write(self, text):
        self.handle.write(text)

    def close(self):
        self.handle.close()

    def delete(self):
        os.unlink(self.name)

def setup_database(con, cur):
    """Create a local database for storing experiment metadata."""
    cur.execute("""
        CREATE TABLE IF NOT EXISTS data (
            step TEXT, 
            path TEXT, 
            parameter TEXT, 
            value TEXT,
            PRIMARY KEY (step, path, parameter))""")
    cur.execute("""
        CREATE TABLE IF NOT EXISTS md5 (
            path TEXT PRIMARY KEY,
            checksum TEXT)""")
    cur.execute("""
        CREATE TABLE IF NOT EXISTS dependencies (
            path TEXT,
            dep_path TEXT,
            dep_checksum TEXT,
            PRIMARY KEY (path, dep_path))
    """)
    con.commit()

def iter_steps(expt):
    """Visit the steps of an experiment in order.
    
    This exploits the fact that every DAG has a topological ordering, which can
    be found by repeatedly choosing and removing a node with no incoming edges.
    """
    dag = {step: expt["Steps"][step]["Depends"][:] for step in expt["Steps"]}

    for i in range(len(dag)):
        # find a step with no prerequisites, and remove it from the DAG
        step = sorted(dag.items(), key = lambda x: len(x[1]))[0][0]
        del dag[step]

        # remove all its outgoing edges
        for k, v in dag.items():
            try:
                dag[k].remove(step)
            except ValueError:
                pass

        yield step

def validate(expt):
    """Sanitize the experiment specification before running anything."""

    # set a global number of processes and threads, serial to be safe
    global_nproc = expt.get("Processes", 1)
    global_nthread = expt.get("Threads", 1)
    shell = os.environ["SHELL"]

    for step_name in expt["Steps"]:
        step = expt["Steps"][step_name]

        # remove steps with no rule
        if step.get("Rule") is None:
            logging.warning("Step {} has no rule and will not be run".format(step_name))
            del expt["Steps"][step_name]
            continue

        # make sure each rule ends with a newline
        if not step["Rule"].endswith(os.linesep):
            step["Rule"] = "{}\n".format(step["Rule"])

        # set the default interpreter to the SHELL environment variable
        if step.get("Interpreter") is None:
            logging.warning("Step {} has no interpreter, assuming it's {}".format(step_name, shell))
            step["Interpreter"] = shell

        # if there's no Startup instruction, set it to blank
        try:
            step["Startup"] = "{}\n".format(step["Startup"].rstrip())
        except KeyError:
            step["Startup"] = ""

        # if there's no extension, it's going to be .dat
        step["Extension"] = step.get("Extension", "dat")

        # every step should have a number of processes and threads
        step["Processes"] = step.get("Processes", global_nproc)
        step["Threads"] = step.get("Threads", global_nthread)

        # parse walltime (defaults to 5 minutes per task)
        walltime = step.get("Walltime", "00:05:00")
        h, m, s = [int(x) for x in walltime.split(":")]
        step["Walltime"] = datetime.timedelta(hours=h, minutes=m, seconds=s)

        # if there are no parameters and/or exclusions, make them an empty dict
        step["Parameters"] = step.get("Parameters", {})
        step["Exclusions"] = step.get("Exclusions", {})

        # make sure every step has prerequisites, and that these are a list (they
        # may be passed as a space-delimited string)
        try:
            step["Depends"] = step["Depends"].split(" ")
        except KeyError:
            step["Depends"] = []
        except AttributeError:
            # prerequisite steps are already a list
            if isinstance(step["Depends"], list):
                pass
            else:
                logging.error("Invalid Depends for step {}".format(step_name))
                logging.error("Depends must be a list or a space-separated string")
                return None

        # make sure the modifications to the steps are saved
        expt["Steps"][step_name] = step

    return expt

def mkdir_p(path):
    """Create a directory, or do nothing error if it already exists."""
    try:
        os.makedirs(path)
    except os.error:
        pass

def iter_parameters(param_dict):
    """Iterate over all possible combinations of parameters.

    Each combination is returned as a dictionary of the same form as the passed
    parameters. It's easiest to illustrate with an example. If we pass in

    param_dict = {
        n: [10, 20],
        p: [0.5, 1.0],
        x: 0
    },

    we will get four iterations, as follows.

    {n: 10, p: 0.5, x: 0}
    {n: 10, p: 1.0, x: 0}
    {n: 20, p: 0.5, x: 0}
    {n: 20, p: 1.0, x: 0}
    """
    for k, v in param_dict.items():
        if isinstance(v, str):
            try:
                v = eval(v)
            except (NameError, SyntaxError):
                pass
        if isinstance(v, (str, int, float)):
            param_dict[k] = [v]
        else:
            param_dict[k] = v

    for param_values in itertools.product(*param_dict.values()):
        yield dict(zip(param_dict.keys(), param_values))

def get_dependencies(expt, step_name, parameters, cur):
    """Find all the prerequisite files for a set of target parameters."""
    step = expt["Steps"][step_name]
    depends = {}

    for prev_step_name in step["Depends"]:
        prev_step = expt["Steps"][prev_step_name]
        common_parameters = set(prev_step["Parameters"].keys()) & set(parameters.keys())

        # if the current step has no parameters in common with the previous step,
        # assume it depends on every file
        if (len(common_parameters) == 0):
            cur.execute("SELECT DISTINCT path FROM data WHERE step = ?", prev_step_name)
            depends[prev_step_name] = list(set(row[0] for row in cur.fetchall()))

        # if there are common parameters, take the subset of files with the same
        # values
        else:
            for i, p in enumerate(common_parameters):
                cur.execute("""
                    SELECT path FROM data WHERE step = ? 
                    AND parameter = ? AND value = ?
                    """, (prev_step_name, p, parameters[p]))
                if i == 0:
                    step_depends = set(row[0] for row in cur.fetchall())
                else:
                    step_depends &= set(row[0] for row in cur.fetchall())
            depends[prev_step_name] = list(step_depends)

    return depends

def needs_update(target, depends, cur):
    """Check if a target needs to be remade."""

    # if the prerequisites are missing, we can't remake the file
    # this isn't as bad as it sounds, because it's normal when there are
    # upstream exclusions
    # so we don't put a warning, just an info
    for step in depends:
        if depends[step] == "":
            logging.info("Prerequisites of step {} for file {} are missing".format(step_name, target))
            return False

    # if the file isn't there at all, we must remake it
    if not os.path.exists(target):
        logging.info("File {} doesn't exist, remaking".format(target))
        return True

    # check if the checksum used to make the target is the same as the current checksum
    for step in depends:
        for file_name in depends[step]:
            cur.execute("SELECT dep_checksum FROM dependencies WHERE path = ? AND dep_path = ?",
                        (target, file_name))
            dep_checksum = cur.fetchone()[0]
            cur.execute("SELECT checksum FROM md5 WHERE path = ?", (file_name,))
            checksum = cur.fetchone()[0]

            if checksum != dep_checksum:
                logging.info("Checksum on dependency {} has changed, remaking file {}".format(file_name, target))
                return True

    logging.info("File {} doesn't need remaking".format(target))
    return False

def find_file(step_dir, step_name, extension, parameters, con, cur):
    """Find a file matching a dictionary of parameters.
    
    If no file exists with those parameters, make a new name for the file by
    hashing the parameters.
    """

    # query the database for a file with all the same parameters

    query = "SELECT path FROM data WHERE step = ? AND parameter = ? AND value = ?"
    query = " INTERSECT ".join([query] * len(parameters))
    query_params = []
    for k, v in parameters.items():
        query_params.extend([step_name, k, v])
    cur.execute(query, tuple(query_params))

    # if there is a matching file in the database, use that
    result = cur.fetchone()
    if result is not None:
        return result[0]

    # if there is no matching file, create a new name and insert it into the
    # database
    basename = hashlib.md5(bytes(str(parameters), "UTF-8")).hexdigest()
    path = os.path.join(step_dir, "{}.{}".format(basename, extension))

    for k, v in parameters.items():
        cur.execute("INSERT INTO data (step, path, parameter, value) VALUES (?, ?, ?, ?)",
                    (step_name, path, k, v))
    con.commit()
    return path

def is_excluded(parameters, exclusions):
    """Check whether a particular set of parameters should be excluded.
    
    This isn't very intelligent, it just checks every exclusion one-by-one.
    """
    param_set = set(parameters.items())
    for excl in exclusions:
        if len(param_set & excl.items()) == len(excl):
            return True
    return False

def store_dependency_checksums(step_name, path, depends, con, cur):
    """Store the checksums of files used to make a target in the database"""

    for dep_step in depends:
        for dep_path in depends[dep_step]:
            cur.execute("SELECT checksum FROM md5 WHERE path = ?", (dep_path,))
            dep_checksum = cur.fetchone()[0]
            cur.execute("""
                INSERT OR REPLACE INTO dependencies (path, dep_path, dep_checksum)
                VALUES (?, ?, ?)
            """, (path, dep_path, dep_checksum))
    con.commit()

def store_checksum(path, con, cur):
    """Store the MD5 checksum for a file in the database."""
    with open(path, "rb") as f:
        checksum = hashlib.md5(f.read()).hexdigest()
        cur.execute("INSERT OR REPLACE INTO md5 (path, checksum) VALUES (?, ?)", 
                    (path, checksum))
        con.commit()

def run_qsub(drivers, targets, con, cur):
    """Run a bunch of driver scripts via qsub, and wait for them to finish"""

    # keep track of the job ID assigned to each driver script
    ids = []
    for driver in drivers:
        proc = subprocess.Popen(["qsub", "-V", driver.name], stdout=subprocess.PIPE)
        ids.append(str(proc.communicate()[0], "UTF-8").split(".")[0])

    # "done" keeps track of which jobs have finished
    done = [False] * len(ids)

    # "target_cur" keeps track of which target is currently being made by each
    # script
    target_cur = [0] * len(ids)

    while not all(done):
        # use qstat to check for running jobs
        proc = subprocess.Popen(["qstat", "-u", os.environ["USER"]], stdout=subprocess.PIPE)
        out = str(proc.communicate()[0], "UTF-8").split("\n")
        out = [line.split() for line in out if len(line.split()) > 9]

        # if a job ID shows up in the qstat output with a status other than
        # "C", the job isn't done
        done = [True] * len(ids)
        for i, id in enumerate(ids):
            for line in out:
                if line[9] != "C" and id in line[0]:
                    done[i] = False
                    break

        # now go through the targets currently being made
        for i, target_list in enumerate(targets):

            # if the job is finished, all its targets are also finished
            if done[i]:
                for j in range(target_cur[i], len(targets[i])):
                    store_checksum(target_list[j], con, cur)
                    target_cur[i] += 1

            # if we're not currently making the last target for this file, and
            # the next target exists, that means the current target is ready
            elif target_cur[i] < len(target_list) - 1 and os.path.exists(target_list[target_cur[i]+1]):
                store_checksum(target_list[target_cur[i]], con, cur)
                target_cur[i] += 1
                
        if not all(done):
            time.sleep(60)

def run_shell(drivers, targets, con, cur):
    """Run a bunch of driver scripts directly in the shell"""
    procs = []
    for driver in drivers:
        procs.append(subprocess.Popen([os.environ["SHELL"], driver.name]))

    # "done" keeps track of which jobs have finished
    done = [False] * len(procs)

    # "target_cur" keeps track of which target is currently being made by each
    # script
    target_cur = [0] * len(procs)

    while not all(done):

        # check if the jobs are finished yet
        for i, proc in enumerate(procs):
            done[i] = proc.poll() is not None

        # now go through the targets currently being made
        for i, target_list in enumerate(targets):

            # if the job is finished, all its targets are also finished
            if done[i]:
                for j in range(target_cur[i], len(targets[i])):
                    store_checksum(target_list[j], con, cur)
                    target_cur[i] += 1

            # if we're not currently making the last target for this file, and
            # the next target exists, that means the current target is ready
            elif target_cur[i] < len(target_list) - 1 and os.path.exists(target_list[target_cur[i]+1]):
                store_checksum(target_list[target_cur[i]], con, cur)
                target_cur[i] += 1
                
        if not all(done):
            time.sleep(60)

def run_step(expt, step_name, con, cur, tmpdir, qsub):
    """Run an experiment step."""

    expt_dir = expt["Name"]
    step_dir = os.path.join(expt_dir, step_name)
    mkdir_p(step_dir)

    step = expt["Steps"][step_name]

    # "files" is a list of which file should be made by each process
    files = []

    # "scripts" is a list of scripts containing the rules to create the files
    scripts = []

    # get the list of excluded parameters
    exclusions = []
    for excl in step["Exclusions"]:
        exclusions.extend(iter_parameters(excl))

    proc_cur = 0
    for parameters in iter_parameters(step["Parameters"]):
        if is_excluded(parameters, exclusions): continue

        # find the file we're making and its dependencies
        target = find_file(step_dir, step_name, step["Extension"], parameters, con, cur)
        depends = get_dependencies(expt, step_name, parameters, cur)

        # check if the file needs updating
        if needs_update(target, depends, cur):

            # open a new script if we need to
            if len(scripts) == proc_cur:
                scripts.append(Tempfile(tmpdir, "w"))
                scripts[proc_cur].write(step["Startup"])
                files.append([])

            files[proc_cur].append(target)

            # insert checksums of dependencies
            store_dependency_checksums(step_name, target, depends, con, cur)

            # add some extras for formatting the rules, plus the names of the
            # target and dependencies
            parameters["yaml"] = yaml.dump(parameters, width=1000).rstrip()
            parameters["seed"] = random.randrange(2**31)
            parameters["$#"] = sum(len(d) for d in depends.values())
            parameters[step_name] = target
            for dep_step in depends:
                parameters[dep_step] = " ".join(depends[dep_step])

            # write the rule to create the file into a script
            command = step["Rule"].format(**parameters)
            logging.debug(command)
            scripts[proc_cur].write(command)

            proc_cur += 1
            if step["Processes"] > 0:
                proc_cur %= step["Processes"]
            con.commit()

    # close the scripts
    for script in scripts:
        script.close()

    # now make a driver to run each of the scripts
    # this seems redundant but it's necessary if the interpreter isn't bash and we're
    # running things with qsub
    drivers = [Tempfile(tmpdir, "w") for s in enumerate(scripts)]
    for i, f in enumerate(drivers):

        f.write("#!{}\n".format(os.environ["SHELL"]))

        # PBS directives
        if qsub:
            walltime = step["Walltime"] * len(files[i])
            hours = int(walltime.seconds / 3600)
            minutes = int((walltime.seconds % 3600) / 60)
            seconds = walltime.seconds % 60

            f.write("#PBS -S {}\n".format(os.environ["SHELL"]))
            f.write("#PBS -l nodes=1:ppn={}\n".format(step["Threads"]))
            f.write("#PBS -l walltime={}:{}:{}\n\n".format(hours, minutes, seconds))
            f.write("cd $PBS_O_WORKDIR\n")

        f.write("{} < {}\n".format(step["Interpreter"], scripts[i].name))
        f.close()

    # submit them all! yay!
    if qsub:
        run_qsub(drivers, files, con, cur)
    else:
        run_shell(drivers, files, con, cur)

    # delete all the scripts
    for f in scripts + drivers:
        f.delete()

def sanitize_database(expt_dir, expt, con, cur):
    """Clean the database."""

    # make sure all the tables are in sync
    cur.execute("SELECT path FROM md5")
    md5_paths = set([row[0] for row in cur.fetchall()])

    cur.execute("SELECT DISTINCT path FROM data")
    data_paths = set([row[0] for row in cur.fetchall()])

    cur.execute("SELECT DISTINCT path FROM dependencies")
    dep_paths = set([row[0] for row in cur.fetchall()])

    for path in md5_paths | data_paths | dep_paths:
        if not os.path.exists(path):
            logging.debug("Deleting file {} from the database (not on filesystem)".format(path))
            cur.execute("DELETE from md5 WHERE path = ?", (path,))
            cur.execute("DELETE from data WHERE path = ?", (path,))
            cur.execute("DELETE from dependencies WHERE path = ?", (path,))
        else:
            with open(path, "rb") as f:
                checksum = hashlib.md5(f.read()).hexdigest()
            cur.execute("SELECT checksum FROM md5 WHERE path = ?", (path,))
            if checksum != cur.fetchone()[0]:
                logging.debug("Deleting file {} from the database (checksum doesn't match)".format(path))
                cur.execute("DELETE from md5 WHERE path = ?", (path,))
                cur.execute("DELETE from data WHERE path = ?", (path,))
                cur.execute("DELETE from dependencies WHERE path = ?", (path,))

    for path in (md5_paths ^ data_paths) | (dep_paths - md5_paths - data_paths):
        logging.debug("Deleting file {} from the database (tables out of sync)".format(path))
        cur.execute("DELETE from md5 WHERE path = ?", (path,))
        cur.execute("DELETE from data WHERE path = ?", (path,))
        cur.execute("DELETE from dependencies WHERE path = ?", (path,))

    # remove files with wrong sets of parameters
    for step_name in expt["Steps"]:

        # get all combinations of parameters which will be used with this step
        params = list(iter_parameters(expt["Steps"][step_name]["Parameters"]))
        exclusions = []
        for excl in expt["Steps"][step_name]["Exclusions"]:
            exclusions.extend(iter_parameters(excl))

        keeps = []
        for p in params:
            if not is_excluded(p, exclusions):
                keeps.append(p)
        params = keeps

        for i, p in enumerate(params):
            params[i] = dict([(k, str(p[k])) for k in sorted(p)])
        params = [set(p.items()) for p in params]

        # get all paths in the database from this step
        cur.execute("SELECT DISTINCT path FROM data WHERE step = ?", (step_name,))
        paths = set([row[0] for row in cur.fetchall()])

        # for each path, see if there is a matching set of parameters
        # if not, delete the path
        for path in paths:
            cur.execute("SELECT parameter, value FROM data WHERE path = ?", (path,))
            cur_params = set(cur.fetchall())

            match = False
            for p in params:
                if len(p) == len(cur_params) and len(cur_params & p) == len(p):
                    match = True
                    break
            if not match:
                logging.debug("Deleting file {} from the database (unused parameters)".format(path))
                cur.execute("DELETE from md5 WHERE path = ?", (path,))
                cur.execute("DELETE from data WHERE path = ?", (path,))
                cur.execute("DELETE from dependencies WHERE path = ?", (path,))

    con.commit()

def sanitize_filesystem(expt_dir, expt, con, cur):
    """Clean the file system."""
    for dirpath, dirnames, filenames in os.walk(expt_dir):
        for basename in filenames:
            path = os.path.join(dirpath, basename)
            cur.execute("SELECT * FROM md5 WHERE path = ?", (path,))
            if cur.fetchone() is None:
                logging.debug("Deleting file {} from filesystem (not in database)".format(path))
                os.unlink(path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a computational experiment.")
    parser.add_argument("yaml_file", type=argparse.FileType("r"))
    parser.add_argument("-q", "--qsub", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-s", "--step")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    expt = validate(yaml.load(args.yaml_file))
    if expt is None:
        sys.exit(1)

    con = sqlite3.connect("{}.sqlite".format(expt["Name"]))
    cur = con.cursor()
    setup_database(con, cur)
    
    expt_dir = expt["Name"]
    mkdir_p(expt_dir)

    tmpdir = ".expttmp"
    mkdir_p(tmpdir)

    sanitize_database(expt_dir, expt, con, cur)
    sanitize_filesystem(expt_dir, expt, con, cur)

    for step_name in iter_steps(expt):
        if args.step is not None and args.step != step_name:
            continue
    
        step = expt["Steps"][step_name]
        nproc = step["Processes"]
        logging.info("Starting step {}".format(step_name))
        run_step(expt, step_name, con, cur, tmpdir, args.qsub)
        con.commit()

    sanitize_database(expt_dir, expt, con, cur)
    sanitize_filesystem(expt_dir, expt, con, cur)
    con.commit()
    con.close()
