#!/bin/sh
# use this script (which requires patchelf utility) to
# fix the paths of the Feff8 shared libraries and binaries

install_dir=$1

base_release='CentOS Linux release 7'

if grep -q "$base_release" /etc/system-release ; then
    echo "will fix binaries"
else
    echo "will only fix binaries on $base_release"
    exit 1
fi

#
# copy shared libs from gcc 4.9.2
cp -pr /usr/local/lib64/libgcc_s.so.1         ./libgcc_s_f8.so
cp -pr /usr/local/lib64/libquadmath.so.0.0.0  ./libquadmath_f8.so
cp -pr /usr/local/lib64/libgfortran.so.3.0.0  ./libgfortran_f8.so

if test -f ./libgfortran_f8.so ; then
    echo "copied gfortran share libaries"
else
    echo "could not copy libgcc_s, libquadmath, and libgfortran"
    exit 1
fi

#
# copy binaries and shared libraries from install location
cp -pr $install_dir/bin/* .
cp -pr $install_dir/lib/* .
mv feff8l.sh feff8l

if test -f ./feff6l ; then
    echo "copied feff binaries and shared libaries"
else
    echo "could not copy installed binaries and shared libraries"
    exit 1
fi
chmod +x lib*.so feff*

#
# fix internal names
patchelf --set-soname libgcc_s_f8.so libgcc_s_f8.so
patchelf --set-soname libgfortran_f8.so libgfortran_f8.so
patchelf --set-soname libquadmath_f8.so libquadmath_f8.so

# tell shared libs to look in '.' for other shared libs
patchelf --set-rpath '$ORIGIN' libgcc_s_f8.so
patchelf --set-rpath '$ORIGIN' libquadmath_f8.so
patchelf --set-rpath '$ORIGIN' libgfortran_f8.so

# fix names to use local copies of shared libs
patchelf --replace-needed libgcc_s.so.1 libgcc_s_f8.so libgfortran_f8.so
patchelf --replace-needed libquadmath.so.0 libquadmath_f8.so libgfortran_f8.so

echo 'fixing Feff8 shared libraries'
for F in libpotph.so libonepath.so libfeff8lpotph.so libfeff8lpath.so; do
    patchelf --set-rpath '$ORIGIN' $F
    patchelf --replace-needed libgcc_s.so.1 libgcc_s_f8.so $F
    patchelf --replace-needed libquadmath.so.0 libquadmath_f8.so $F
    patchelf --replace-needed libgfortran.so.3 libgfortran_f8.so $F
done

echo 'fixing Feff8 binaries'
for F in (feff6l feff8l_ff2x feff8l_genfmt feff8l_pathfinder
	  feff8l_pot feff8l_rdinp feff8l_xsph); do
  patchelf --set-rpath '$ORIGIN' $F
  patchelf --replace-needed libgcc_s.so.1 libgcc_s_f8.so $F
  patchelf --replace-needed libquadmath.so.0 libquadmath_f8.so $F
  patchelf --replace-needed libgfortran.so.3 libgfortran_f8.so $F

done
