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

try:
    subprocess_run = subprocess.run
except AttributeError:
    subprocess_run = subprocess.call

def mkdir_p(path):
    try:
        os.makedirs(path)
    except os.error:
        pass

def iter_parameters(param_dict):
    for k, v in param_dict.items():
        if isinstance(v, (str, int, float)):
            param_dict[k] = [v]

    for param_values in itertools.product(*param_dict.values()):
        yield dict(zip(param_dict.keys(), param_values))

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
        common_parameters = set(prev_step["Parameters"].keys()) & set(parameters.keys())

        # if the current step has no parameters in common with the previous step,
        # assume it depends on every file
        if (len(common_parameters) == 0):
            cur.execute("SELECT file FROM data WHERE experiment = ? AND step = ?",
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
    if not os.path.exists(target):
        return True
    for step in depends:
        for file_name in depends[step].split(" "):
            with open(file_name, "rb") as f:
                checksum = hashlib.md5(f.read()).hexdigest()
                cur.execute("SELECT checksum FROM md5 WHERE file = ?", (file_name,))
                row = cur.fetchone()
                if row is None or row[0] != checksum:
                    return True
    return False

def run_step(expt, step_name, con, cur, nproc):
    step = expt["Steps"][step_name]

    interpreter = shlex.split(step["Interpreter"])
    try:
        startup = step["Startup"].rstrip() + "\n"
    except KeyError:
        startup = ""
    startup = bytes(startup, "UTF-8")

    # steps are placed in expt_name/step_name
    expt_name = expt["Name"]
    step_dir = os.path.join(expt_name, step_name)
    mkdir_p(step_dir)

    # keep track of all the files we made
    files = []

    # make a file for each combination of parameters
    proc = [0] * nproc
    proc_cur = 0

    for parameters in iter_parameters(step["Parameters"]):

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

        if len(old_basename) == 0:
            basename = hashlib.md5(bytes(str(parameters), "UTF-8")).hexdigest()[:8]
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

            # start the interpreter if it's not already open
            if proc[proc_cur] == 0:
                proc[proc_cur] = subprocess.Popen(interpreter, stdin=subprocess.PIPE)
                proc[proc_cur].stdin.write(startup)

            # insert file and parameters into database
            cur.executemany("""
                INSERT OR REPLACE INTO data (experiment, step, file, parameter, value)
                VALUES (?, ?, ?, ?, ?)
            """, ((expt_name, step_name, basename, k, v) for k, v in parameters.items()))

            # add some extras for formatting the rules
            parameters["yaml"] = yaml.dump(parameters, width=1000).rstrip()
            parameters[step_name] = target
            parameters.update(depends)
            parameters["$#"] = sum(len(d.split(" ")) for d in depends.values())

            # create the file
            command = step["Rule"].format(**parameters)
            if not command.endswith(os.linesep):
                command = "{}\n".format(command)
            proc[proc_cur].stdin.write(bytes(command, "UTF-8"))
            files.append(target)

            proc_cur = (proc_cur + 1) % nproc

            # commit changes to the database
            con.commit()

    for i in range(nproc):
        try:
            proc[i].communicate()
        except AttributeError:
            pass

    for target in files:
        with open(target, "rb") as f:
            checksum = hashlib.md5(f.read()).hexdigest()
            cur.execute("INSERT OR REPLACE INTO md5 (file, checksum) VALUES (?, ?)", (target, checksum))
    con.commit()

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
            file TEXT PRIMARY KEY,
            checksum TEXT)""")

def iter_steps(expt):
    dag = {step: expt["Steps"][step].get("Depends") for step in expt["Steps"]}
    for k, v in dag.items():
        if not v:
            dag[k] = []
        elif isinstance(v, str):
            dag[k] = [v]

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
    yaml_file = sys.argv[1]

    with open(yaml_file) as f:
        expt = yaml.load(f)
    con = sqlite3.connect("expt.sqlite")
    cur = con.cursor()
    setup_database(con, cur)
    
    expt_dir = expt["Name"]
    mkdir_p(expt_dir)

    try:
        nproc = expt["Processes"]
    except KeyError:
        nproc = 1
    
    for step_name in iter_steps(expt):
        try:
            expt["Steps"][step_name]["Depends"] = expt["Steps"][step_name]["Depends"].split(" ")
        except AttributeError:
            pass
        run_step(expt, step_name, con, cur, nproc)

    con.close()
