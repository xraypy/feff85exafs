#!/usr/bin/python

import sys
import subprocess
import glob
import re
from   os      import unlink, makedirs, chdir
from   os.path import isdir,  join
from   shutil  import rmtree
import pystache
import json

## write the feff.inp file
target = join(sys.argv[1], 'baseline')
if isdir(target):
    rmtree(target)
makedirs(target)
renderer = pystache.Renderer()

with open(join(target,'feff.inp'), 'w') as inp:
    inp.write(renderer.render_path( join(sys.argv[1], sys.argv[1] + '.mustache'),                # material/material.mustache
                                    json.load(open(join(sys.argv[1], sys.argv[1] + '.json'))) )) # material/material.json

chdir(target)

## run the f85e shell script, which emulates the behavior of the monolithic Feff application
subprocess.call('/home/bruce/git/feff85exafs/bin/f85e');

## cull the files we don't need for testing
feffoutput = glob.glob("*")
for f in sorted(feffoutput):
    tosave = re.compile("feff(\d+\.dat|\.inp)|(chi|files|paths|xmu)\.dat|f85e.log")
    if not tosave.match(f):
        unlink(f)
