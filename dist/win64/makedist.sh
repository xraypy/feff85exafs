#!/bin/sh
# use this script with MSYS2 / MinGW to copy needed dlls

install_dir=$1

cp $install_dir/bin/* .
cp $install_dir/lib/* .

cp /mingw64/bin/libgcc_s_seh-1.dll .
cp /mingw64/bin/libgfortran-5.dll .
cp /mingw64/bin/libquadmath-0.dll .
cp /mingw64/bin/libwinpthread-1.dll .

