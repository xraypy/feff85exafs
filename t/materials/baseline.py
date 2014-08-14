#!/usr/bin/python

import sys, subprocess, glob, re
from   os      import unlink, makedirs, chdir
from   os.path import isdir,  realpath, join
from   shutil  import rmtree
import pystache, json


repotop = realpath(join('..','..'))
folder = sys.argv[1]
if folder[-1] == '/': folder = folder[:-1]

## write the feff.inp file
target = join(folder, 'baseline')
if isdir(target): rmtree(target)
makedirs(target)
renderer = pystache.Renderer()

with open(join(target,'feff.inp'), 'w') as inp:
    inp.write(renderer.render_path( join(folder, folder + '.mustache'),                # material/material.mustache
                                    json.load(open(join(folder, folder + '.json'))) )) # material/material.json

chdir(target)

f85escript = join(repotop, 'bin', 'f85e')
## run the f85e shell script, which emulates the behavior of the monolithic Feff application
subprocess.call(f85escript);

## cull the files we don't need for testing
feffoutput = glob.glob("*")
for f in sorted(feffoutput):
    tosave = re.compile("feff(\d+\.dat|\.inp)|(chi|files|paths|xmu)\.dat|f85e.log")
    if not tosave.match(f): unlink(f)
