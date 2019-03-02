#!/usr/bin/env python 

from __future__ import print_function

import argparse
import pandas as pd
import yaml
import json


parser = argparse.ArgumentParser()

parser.add_argument("--run", action="store", nargs="?")
parser.add_argument("fname")

# validation
args = parser.parse_args()
fname = args.fname
# print(args)  # debug

df = pd.read_csv(fname)

if args.run is not None:
    run = args.run

    df.run = df.run.astype(str)
    constants = df[df.run == run].squeeze().to_json(double_precision=15)
    constants = json.loads(constants)

    print(yaml.dump(constants, default_flow_style=False))


