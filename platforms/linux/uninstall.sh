#!/bin/sh

files='/usr/local/share/larch/dlls/linux32/_feffpathwrapper.so
/usr/local/share/larch/plugins/xafs/scatteringpath.py
/usr/local/share/larch/modules/feffpathwrapper.py
/usr/local/bin/ff2x
/usr/local/bin/genfmt
/usr/local/bin/pathfinder
/usr/local/bin/pot
/usr/local/bin/rdinp
/usr/local/bin/xsph
/usr/local/lib/libfeffpath.so
/usr/local/lib/libonepath.so
/usr/local/lib/json_module.mod
/usr/local/lib/libjsonfortran.a
'

for f in $files; do
    rm -fv $f
done
