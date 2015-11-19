#!/usr/bin/env python3

import shlex
import subprocess
import itertools
import yaml
import os
import hashlib
import sys
from pprint import pprint
import random
import sqlite3
import hashlib
import logging
import random
import tempfile
import time

TMPDIR = ".expttmp"

try:
    subprocess_run = subprocess.run
except AttributeError:
    subprocess_run = subprocess.call

def mkdir_p(path):
    try:
        os.makedirs(path)
    except os.error:
        pass

def iter_parameters(param_dict, exclusions):
    for k, v in param_dict.items():
        if isinstance(v, str):
            try:
                param_dict[k] = eval(v)
                v = param_dict[k]
            except:
                pass
        if isinstance(v, (str, int, float)):
            param_dict[k] = [v]

    for param_values in itertools.product(*param_dict.values()):
        pdict = dict(zip(param_dict.keys(), param_values))
        keep = True
        for excl in exclusions:
            shared_items = set(pdict.items()) & set(excl.items())
            if len(shared_items) == len(excl.items()):
                keep = False
        if keep:
            yield pdict

def get_dependencies(expt, step_name, parameters, cur):
    step = expt["Steps"][step_name]
    depends = {}

    try:
        if not step["Depends"]:
            return depends
        elif isinstance(step["Depends"], str):
            prev_steps = [step["Depends"]]
        else:
            prev_steps = step["Depends"]
    except KeyError:
        return depends

    for prev_step_name in prev_steps:
        prev_step = expt["Steps"][prev_step_name]
        try:
            common_parameters = set(prev_step["Parameters"].keys()) & set(parameters.keys())
        except KeyError:
            common_parameters = []

        # if the current step has no parameters in common with the previous step,
        # assume it depends on every file
        if (len(common_parameters) == 0):
            cur.execute("SELECT file FROM md5 WHERE experiment = ? AND step = ?",
                        (expt["Name"], prev_step_name))
            depends[prev_step_name] = list(set(row[0] for row in cur.fetchall()))

        # if there are common parameters, take the subset of files with the same
        # values
        else:
            for i, p in enumerate(common_parameters):
                cur.execute("""
                    SELECT file FROM data WHERE experiment = ? 
                    AND step = ? AND parameter = ? AND value = ?
                    """, (expt["Name"], prev_step_name, p, parameters[p]))
                if i == 0:
                    step_depends = set(row[0] for row in cur.fetchall())
                else:
                    step_depends &= set(row[0] for row in cur.fetchall())
            depends[prev_step_name] = list(step_depends)

        extension = prev_step["Extension"]
        path = os.path.join(expt["Name"], prev_step_name)

        depends[prev_step_name] = sorted(depends[prev_step_name])
        for i, d in enumerate(depends[prev_step_name]):
            depends[prev_step_name][i] = os.path.join(path, "{}.{}".format(d, extension))
        depends[prev_step_name] = " ".join(depends[prev_step_name])
    return depends

def needs_update(target, depends, cur):
    for step in depends:
        if depends[step] == "":
            logging.info("Prerequisites of step {} for file {} are missing".format(step_name, target))
            return False

    if not os.path.exists(target):
        logging.info("File {} doesn't exist, remaking".format(target))
        return True

    for step in depends:
        for file_name in depends[step].split(" "):
            with open(file_name, "rb") as f:
                checksum = hashlib.md5(f.read()).hexdigest()
                cur.execute("SELECT checksum FROM md5 WHERE path = ?", (file_name,))
                row = cur.fetchone()
                if row is None:
                    logging.error("Dependency {} doesn't exist".format(file_name))
                    return False
                if row[0] != checksum:
                    logging.info("Checksum on dependency {} has changed, remaking file {}".format(file_name, target))
                    return True

    logging.info("File {} doesn't need remaking".format(target))
    return False

def run_step(expt, step_name, con, cur, nproc):
    step = expt["Steps"][step_name]

    # if there's no interpreter, nothing to do
    try:
        interpreter = step["Interpreter"]
    except KeyError:
        return True

    try:
        startup = step["Startup"].rstrip() + "\n"
    except KeyError:
        startup = ""
    if sys.version_info.major == 3:
        startup = bytes(startup, "UTF-8")

    # steps are placed in expt_name/step_name
    expt_name = expt["Name"]
    step_dir = os.path.join(expt_name, step_name)
    mkdir_p(step_dir)

    # keep track of all the files we made
    files = []

    scripts = [0] * nproc
    ntasks = [0] * nproc
    proc_cur = 0

    exclusions = []
    try:
        for excl in step["Exclusions"]:
            exclusions.extend(list(iter_parameters(excl, [])))
    except KeyError:
        pass

    for parameters in iter_parameters(step.get("Parameters"), exclusions):

        # check if a file with these parameters is there already
        old_basename = None
        for k, v in parameters.items():
            cur.execute("""
                SELECT file FROM data WHERE experiment = ? AND step = ? 
                AND parameter = ? AND value = ?
                """, (expt_name, step_name, k, v))
            if old_basename is None:
                old_basename = set([row[0] for row in cur.fetchall()])
            else:
                old_basename &= set([row[0] for row in cur.fetchall()])

        # no parameters, should be an external file
        if old_basename is None:
            try:
                old_basename = [step["ExternalFile"]]
            except KeyError:
                raise RuntimeError("No parameters and no external file for step {}\n".format(step_name))

        if len(old_basename) == 0:
            if sys.version_info.major == 3:
                basename = hashlib.md5(bytes(str(parameters), "UTF-8")).hexdigest()[:8]
            else:
                basename = hashlib.md5(str(parameters)).hexdigest()[:8]
        elif len(old_basename) == 1:
            basename = old_basename.pop()
        else:
            raise RuntimeError("Database integrity lost: redundant files "
            "for step {} of experiment {}\n".format(step_name, expt_name))

        file_name = "{}.{}".format(basename, step["Extension"])
        target = os.path.join(step_dir, file_name)

        # gather dependencies
        depends = get_dependencies(expt, step_name, parameters, cur)

        # check if the file needs updating
        if (needs_update(target, depends, cur)):
            if nproc <= 0:
                tf = tempfile.mkstemp(dir=TMPDIR)
                scripts.append((os.fdopen(tf[0], "w"), tf[1]))
                scripts[-1][0].write(str(startup, "UTF-8"))
                
            elif scripts[proc_cur] == 0:
                tf = tempfile.mkstemp(dir=TMPDIR)
                scripts[proc_cur] = (os.fdopen(tf[0], "w"), tf[1])
                scripts[proc_cur][0].write(str(startup, "UTF-8"))

            # insert file and parameters into database
            cur.executemany("""
                INSERT OR REPLACE INTO data (experiment, step, file, parameter, value)
                VALUES (?, ?, ?, ?, ?)
            """, ((expt_name, step_name, basename, k, v) for k, v in parameters.items()))

            # add some extras for formatting the rules
            parameters["seed"] = random.randrange(2**31)
            parameters["yaml"] = yaml.dump(parameters, width=1000).rstrip()
            parameters[step_name] = target
            parameters.update(depends)
            parameters["$#"] = sum(len(d.split(" ")) for d in depends.values())

            # create the file
            command = step["Rule"].format(**parameters)
            if not command.endswith(os.linesep):
                command = "{}\n".format(command)
            logging.debug(command)
            
            if nproc <= 0:
                scripts[-1][0].write(command)
            else:
                scripts[proc_cur][0].write(command)
                ntasks[proc_cur] += 1
            files.append((basename, target))

            proc_cur = (proc_cur + 1) % nproc

            # commit changes to the database
            con.commit()

    jobscripts = [0] * len(scripts)
    ids = []
    for i in range(len(scripts)):
        if scripts[i] != 0:
            scripts[i][0].write("\n")
            scripts[i][0].close()
            tf = tempfile.mkstemp(dir=TMPDIR)
            jobscripts[i] = (os.fdopen(tf[0], "w"), tf[1])
            
            jobscripts[i][0].write("#!/bin/bash\n")
            jobscripts[i][0].write("#PBS -S /bin/bash\n")
            try:
                jobscripts[i][0].write("#PBS -l nodes=1:ppn={}\n".format(step["Threads"]))
            except KeyError:
                jobscripts[i][0].write("#PBS -l nodes=1:ppn=1\n")
            try:
                jobscripts[i][0].write("#PBS -l walltime={}:00:00\n".format(step["Walltime"]))
            except KeyError:
                jobscripts[i][0].write("#PBS -l walltime={}:00:00\n".format(ntasks[i]))
            jobscripts[i][0].write("\n")

            jobscripts[i][0].write("cd $PBS_O_WORKDIR\n")
            jobscripts[i][0].write("{} < {}\n".format(interpreter, scripts[i][1]))
            jobscripts[i][0].close()

            proc = subprocess.Popen(["qsub", "-V", jobscripts[i][1]], stdout=subprocess.PIPE)
            ids.append(str(proc.communicate()[0], "UTF-8").split(".")[0])
            
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

    for _, fn in [x for x in scripts + jobscripts if x != 0]:
        os.unlink(fn)

    # add new files to DB
    ok = True
    for basename, target in files:
        if not os.path.exists(target):
            logging.error("File {} was not created".format(target))
            ok = False
            continue
        with open(target, "rb") as f:
            checksum = hashlib.md5(f.read()).hexdigest()
            cur.execute("INSERT OR REPLACE INTO md5 (experiment, step, file, path, checksum) VALUES (?, ?, ?, ?, ?)", (expt_name, step_name, basename, target, checksum))
    con.commit()
    return ok

def setup_database(con, cur):
    cur.execute("""
        CREATE TABLE IF NOT EXISTS data (
            experiment TEXT, 
            step TEXT, 
            file TEXT, 
            parameter TEXT, 
            value TEXT,
            PRIMARY KEY (experiment, step, file, parameter))""")
    cur.execute("""
        CREATE TABLE IF NOT EXISTS md5 (
            experiment TEXT, 
            step TEXT, 
            file TEXT,
            path TEXT PRIMARY KEY,
            checksum TEXT)""")

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
                    cur.execute("INSERT OR REPLACE INTO md5 (experiment, step, file, path, checksum) VALUES (?, ?, ?, ?, ?)", (expt, step, basename, path, checksum))
    con.commit()

    # remove files which don't exist from the DB
    cur.execute("SELECT path FROM md5 WHERE experiment = ?", (expt,))
    for row in cur.fetchall():
        if not os.path.exists(row[0]):
            logging.info("Deleting file {}, which is not on the filesystem, from the DB".format(row[0]))
            cur.execute("DELETE FROM md5 WHERE path = ?", row)
    con.commit()

    # remove files not in the DB from the filesystem
    cur.execute("SELECT DISTINCT experiment, step FROM data WHERE experiment = ?", (expt,))
    for row in cur.fetchall():
        path = os.path.join(*row)
        if not os.path.exists(path):
            continue
        for f in os.listdir(path):
            cur.execute("SELECT checksum FROM md5 WHERE path = ?", (os.path.join(path, f),))
            if cur.fetchone() is None:
                logging.info("Deleting file {}, which is not in the DB, from the filesystem".format(os.path.join(path, f)))
                os.remove(os.path.join(path, f))

    cur.execute("SELECT DISTINCT experiment, step, file FROM data NATURAL LEFT OUTER JOIN md5 WHERE path IS NULL AND experiment = ?", (expt,))
    for row in cur.fetchall():
        logging.info("Deleting file {}/{}/{}, which is not in md5 table, from data table".format(*row))
        cur.execute("DELETE FROM data WHERE experiment = ? AND step = ? AND file = ?", row)

def iter_steps(expt):
    dag = {step: expt["Steps"][step].get("Depends") for step in expt["Steps"]}
    for k, v in dag.items():
        if not v:
            dag[k] = []
        elif isinstance(v, str):
            dag[k] = v.split(" ")

    for i in range(len(dag)):
        step = sorted(dag.items(), key = lambda x: len(x[1]))[0][0]
        del dag[step]
        for k, v in dag.items():
            try:
                dag[k].remove(step)
            except ValueError:
                pass
        yield step

if __name__ == "__main__":
    # TODO: proper argument parsing
    yaml_file = sys.argv[1]
    logging.basicConfig(level=logging.DEBUG)
    mkdir_p(TMPDIR)

    with open(yaml_file) as f:
        expt = yaml.load(f)
    con = sqlite3.connect("{}.sqlite".format(expt["Name"]))
    cur = con.cursor()
    setup_database(con, cur)
    
    expt_dir = expt["Name"]
    mkdir_p(expt_dir)
    
    sanitize(expt["Name"], con, cur)
    for step_name in iter_steps(expt):
    
        try:
            nproc = expt["Steps"][step_name]["Processes"]
        except KeyError:
            try:
                nproc = expt["Processes"]
            except KeyError:
                nproc = 1
        logging.info("Starting step {}".format(step_name))
        try:
            expt["Steps"][step_name]["Depends"] = expt["Steps"][step_name]["Depends"].split(" ")
        except AttributeError:
            pass
        except KeyError:
            pass
        if not run_step(expt, step_name, con, cur, nproc):
            logging.error("Step {} falied".format(step_name))
        else:
            logging.info("Finished step {}".format(step_name))
        sanitize(expt["Name"], con, cur)

    con.close()
