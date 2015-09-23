import itertools
import re
import yaml
import os.path

def clean_rule(rule):
    rule = re.sub("(\s)\s+", "\\1", rule)
    rule = re.sub("^\s", "", rule)
    rule = re.sub("\s$", "", rule)
    rule = re.sub("\n", "\n\t", rule)
    return rule

def create_makefile(params, targets, target_dir):

    # targets is keyed by tuples of parameters, make them sorted
    targets = {tuple(sorted(k)):v for k, v in targets.items()}

    # clean up the rules so they can be directly printed to a Makefile
    for k in targets:
        for i, rule in enumerate(targets[k]):
            target, depends, command = rule
            targets[k][i] = (target, depends, clean_rule(command))

    # if any of the parameters are all numeric, pad with zeroes so they are all
    # the same width
    for k, v in params.items():
        if all(isinstance(x, (int)) for x in v):
            v = [str(x) for x in v]
            width = max(len(x) for x in v)
            params[k] = ["0"*(width - len(x)) + x for x in v]

        elif all(isinstance(x, (float, int)) for x in v):
            s = [str(x) for x in v]
            left_width = max(len(x.split(".")[0]) for x in s)
            right_width = max(len(x.split(".")[1]) if "." in x else 0 for x in s)
            width = left_width + 1 + right_width
            fmt = "{{:0{width}.{right_width}f}}".format(width=width, right_width=right_width)
            params[k] = [fmt.format(x) for x in v]

    # keep track of a list of targets, and also the text for the Makefile
    to_make = []
    makefile = ""

    # go through targets in order of fewest-to-most parameter dependencies
    for num_params in range(len(params)+1, -1, -1):
        for combo in itertools.combinations(params.keys(), num_params):
            key = tuple(sorted(combo))
            try:
                for target, depends, command in targets[key]:
                    for p in itertools.product(*[params[k] for k in key]):
                        d = dict(zip(key, p))
                        d["yaml"] = yaml.dump(d, width=1000)[:-1]
                        d["file_basename"] = os.path.join(target_dir, "_".join(map(str, p)))
                        makefile += "{}: {}\n\t{}\n\n".format(target.format(**d),
                                depends.format(**d), command.format(**d))
                        to_make.append(target.format(**d))
            except KeyError as e:
                pass
    return makefile, to_make

