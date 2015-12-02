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
            step TEXT, 
            path TEXT PRIMARY KEY,
            checksum TEXT)""")
    cur.execute("""
        CREATE TABLE IF NOT EXISTS dependencies (
            step TEXT,
            path TEXT,
            dep_step TEXT,
            dep_path TEXT,
            dep_checksum TEXT,
            PRIMARY KEY (step, path, dep_step, dep_path))
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
            cur.execute("SELECT path FROM md5 WHERE step = ?", prev_step_name)
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
            cur.execute("SELECT checksum FROM md5 WHERE step = ? AND path = ?", (step, file_name))
            checksum = cur.fetchone()[0]

            if checksum != dep_checksum:
                logging.info("Checksum on dependency {} has changed, remaking file {}".format(file_name, target))
                return True

    logging.info("File {} doesn't need remaking".format(target))
    return False

def run_step(expt, step_name, con, cur, tmpdir, qsub):
    """Run an experiment step."""

    expt_dir = expt["Name"]
    step = expt["Steps"][step_name]
    nproc = step["Processes"]

    # steps are placed in expt_name/step_name
    step_dir = os.path.join(expt_dir, step_name)
    mkdir_p(step_dir)

    # keep track of all the files we made
    files = []
    scripts = []
    ntasks = []
    proc_cur = 0

    # get the list of excluded parameters
    exclusions = []
    for excl in step["Exclusions"]:
        exclusions.extend(iter_parameters(excl))

    for parameters in iter_parameters(step["Parameters"]):
        # if the set of parameters completely matches one of the exclusions, skip it
        exclude = False
        for excl in exclusions:
            if len(set(parameters.items()) & set(excl.items())) == len(excl):
                exclude = True
                break
        if exclude:
            continue

        # check if a file with these parameters is there already
        old_file = None
        for k, v in parameters.items():
            cur.execute("""
                SELECT path FROM data WHERE step = ? 
                AND parameter = ? AND value = ?
                """, (step_name, k, v))
            if old_file is None:
                old_file = set([row[0] for row in cur.fetchall()])
            else:
                old_file &= set([row[0] for row in cur.fetchall()])
                if len(old_file) == 0:
                    break

        # if there's no extant file, create a new one from a hash of the parameters
        if len(old_file) == 0:
            basename = hashlib.md5(bytes(str(parameters), "UTF-8")).hexdigest()
            target = os.path.join(step_dir, "{}.{}".format(basename, step["Extension"]))
        
        # if there is an extant file, just use that
        elif len(old_file) == 1:
            target = old_file.pop()

        # if there is more than one extant file, something has gone wrong
        else:
            logging.error("Database integrity lost")
            logging.error("Redundant files for step {}".format(step_name))
            sys.exit(1)

        # gather dependencies
        depends = get_dependencies(expt, step_name, parameters, cur)

        # check if the file needs updating
        if needs_update(target, depends, cur):

            # open a new script if we need to
            if len(scripts) == proc_cur:
                fd, fn = tempfile.mkstemp(dir=tmpdir)
                f = os.fdopen(fd, "w")
                f.write(step["Startup"])
                scripts.append((f, fn))
                ntasks.append(0)

            # insert file and parameters into database
            cur.executemany("""
                INSERT OR REPLACE INTO data (step, path, parameter, value)
                VALUES (?, ?, ?, ?)
            """, ((step_name, target, k, v) for k, v in parameters.items()))

            # insert checksums of dependencies
            for dep_step in depends:
                for dep_path in depends[dep_step]:
                    cur.execute("SELECT checksum FROM md5 WHERE path = ?", (dep_path,))
                    dep_checksum = cur.fetchone()[0]
                    cur.execute("""
                        INSERT OR REPLACE INTO dependencies (step, path, dep_step, dep_path, dep_checksum)
                        VALUES (?, ?, ?, ?, ?)
                    """, (step_name, target, dep_step, dep_path, dep_checksum))

            # add some extras for formatting the rules
            parameters["yaml"] = yaml.dump(parameters, width=1000).rstrip()
            parameters["seed"] = random.randrange(2**31)
            parameters[step_name] = target
            parameters["$#"] = sum(len(d) for d in depends.values())

            # add all the prerequisite files
            for dep_step in depends:
                depends[dep_step] = " ".join(depends[dep_step])
            parameters.update(depends)

            # write the rule to create the file into a script
            command = step["Rule"].format(**parameters)
            if not command.endswith(os.linesep):
                command = "{}\n".format(command)
            logging.debug(command)
            scripts[proc_cur][0].write(command)
            files.append(target)
            ntasks[proc_cur] += 1

            proc_cur += 1
            if nproc > 0:
                proc_cur %= nproc

            # commit changes to the database
            con.commit()

    # we don't need the file handles anymore, just the names
    for i, script in enumerate(scripts):
        script[0].close()
        scripts[i] = script[1]

    # now make a second script to run each of the lists of commands
    # this seems redundant but it's necessary if the interpreter isn't bash and we're
    # running things with qsub
    run_scripts = []
    for i, script in enumerate(scripts):

        fd, fn = tempfile.mkstemp(dir=tmpdir)
        f = os.fdopen(fd, "w")
        f.write("#!{}\n".format(os.environ["SHELL"]))

        # PBS directives
        if qsub:
            walltime = step["Walltime"] * ntasks[i]
            hours = int(walltime.seconds / 3600)
            minutes = int((walltime.seconds % 3600) / 60)
            seconds = walltime.seconds % 60

            f.write("#PBS -S {}\n".format(os.environ["SHELL"]))
            f.write("#PBS -l nodes=1:ppn={}\n".format(step["Threads"]))
            f.write("#PBS -l walltime={}:{}:{}\n\n".format(hours, minutes, seconds))
            f.write("cd $PBS_O_WORKDIR\n")

        f.write("{} < {}\n".format(step["Interpreter"], script))
        f.close()
        run_scripts.append(fn)

    # submit them all! yay!
    ids = []
    for script in run_scripts:
        if qsub:
            proc = subprocess.Popen(["qsub", "-V", script], stdout=subprocess.PIPE)
            ids.append(str(proc.communicate()[0], "UTF-8").split(".")[0])
        else:
            proc = subprocess.Popen([os.environ["SHELL"], script])
            ids.append(proc)

    # if we submitted the jobs to qsub, we have to monitor the queue to figure out
    # when they are finished
    if qsub:
        done = False
        while not done:
            proc = subprocess.Popen(["qstat", "-u", "rmcclosk"], stdout=subprocess.PIPE)
            out = str(proc.communicate()[0], "UTF-8")

            done = True
            for id in ids:
                for line in out.split("\n"):
                        line = line.split()
                        if len(line) > 9 and line[9] != "C" and id in line[0]:
                            done = False
                            time.sleep(60)
                            break
                if not done:
                    break

    # if we are just running shells, wait for them to finish
    else:
        [id.communicate() for id in ids]

    # delete all the scripts
    for f in scripts + run_scripts:
        os.unlink(f)

    # add checksums for the new files to the database
    for target in files:
        if not os.path.exists(target):
            logging.error("File {} was not created".format(target))
            continue
        with open(target, "rb") as f:
            checksum = hashlib.md5(f.read()).hexdigest()
            cur.execute("INSERT OR REPLACE INTO md5 (step, path, checksum) VALUES (?, ?, ?)", (step_name, target, checksum))
    con.commit()

def sanitize(expt, con, cur):
    # add files which were missed before
    for step in os.listdir(expt):
        if not os.path.exists(os.path.join(expt, step)):
            continue
        for file in os.listdir(os.path.join(expt, step)):
            path = os.path.join(expt, step, file)
            cur.execute("SELECT * FROM md5 WHERE path = ?", (path,))
            if cur.fetchone() is None:
                logging.info("Adding checksum for file {} to DB".format(path))
                basename = file.split(".")[0]
                with open(path, "rb") as f:
                    checksum = hashlib.md5(f.read()).hexdigest()
                    cur.execute("INSERT OR REPLACE INTO md5 (step, path, checksum) VALUES (?, ?, ?)", (step, path, checksum))
    con.commit()

    # remove files which don't exist from the DB
    cur.execute("SELECT path FROM md5")
    for row in cur.fetchall():
        if not os.path.exists(row[0]):
            logging.info("Deleting file {}, which is not on the filesystem, from the DB".format(row[0]))
            cur.execute("DELETE FROM md5 WHERE path = ?", row)
    con.commit()

    # remove files not in the DB from the filesystem
    cur.execute("SELECT DISTINCT step FROM data")
    for row in cur.fetchall():
        path = os.path.join(*row)
        if not os.path.exists(path):
            continue
        for f in os.listdir(path):
            cur.execute("SELECT checksum FROM md5 WHERE path = ?", (os.path.join(path, f),))
            if cur.fetchone() is None:
                logging.info("Deleting file {}, which is not in the DB, from the filesystem".format(os.path.join(path, f)))
                os.remove(os.path.join(path, f))

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
    sanitize(expt_dir, con, cur)

    tmpdir = ".expttmp"
    mkdir_p(tmpdir)

    for step_name in iter_steps(expt):
        if args.step is not None and args.step != step_name:
            continue
    
        step = expt["Steps"][step_name]
        nproc = step["Processes"]
        logging.info("Starting step {}".format(step_name))
        run_step(expt, step_name, con, cur, tmpdir, args.qsub)
        sanitize(expt_dir, con, cur)

    con.close()
