#!/bin/sh
# use this script (which requires install_name_tool utility from
# the Xcode commandlne tools) to fix the paths of the Feff8 shared
# libraries and binaries
#
#
#   mv libca.3.16.2.dylib libca.dylib
#   mv libCom.3.16.2.dylib libComPYEPICS.dylib
#   install_name_tool -id libca.dylib  libca.dylib
#   install_name_tool -id libComPYEPICS.dylib  libComPYEPICS.dylib
#   install_name_tool -add-rpath . libca.dylib
#   install_name_tool -change libComPYEPICS.dylib @loader_path/./libComPYEPICS.dylib libca.dylib

#!/bin/sh

install_dir=$1
install_dir='../../local_install'

gfortran_lib=`gfortran -print-file-name=libgfortran.dylib`
gfortran_dir=`dirname $gfortran_lib`

# copy shared libs from gcc 
cp -pr $gfortran_dir/libgcc_s.1.dylib .
cp -pr $gfortran_dir/libgfortran.5.dylib .
cp -pr $gfortran_dir/libquadmath.0.dylib .

if test -f ./libgfortran.5.dylib ; then
    echo "copied gfortran share libaries"
else
    echo "could not copy libgcc_s, libquadmath, and libgfortran"
    exit 1
fi

cp -pr $install_dir/bin/* .
cp -pr $install_dir/lib/* .
mv feff8l.sh feff8l

for fname in libgfortran.5.dylib libquadmath.0.dylib libgcc_s.1.dylib \
				 libfeff8lpath.dylib libpotph.dylib \
				 libfeff8lpotph.dylib libonepath.dylib \
				 feff6l feff8l_ff2x feff8l_pathfinder \
				 feff8l_rdinp feff8l_genfmt feff8l_pot feff8l_xsph ; do
    chmod 755 $fname
    install_name_tool -id @loader_path/./$lname $fname
    install_name_tool -add_rpath . $fname

    for dname in /usr/local/lib/gcc/10 /usr/local/opt/gcc/lib/gcc/10 \
				       /usr/local/Cellar/gcc/10.2.0/lib/gcc/10 ; do
	for xname in libgfortran.5.dylib libquadmath.0.dylib libgcc_s.1.dylib; do
	    install_name_tool -change $dname/$xname @loader_path/./$xname $fname
	done
    done
done

# # 
# if test -f ./feff6l ; then
#     echo "copied feff binaries and shared libaries"
# else
#     echo "could not copy installed binaries and shared libraries"
#     exit 1
# fi
# chmod +x lib*.so feff*
# 
# echo 'ready ...'


# 
# #
# # fix internal names
# patchelf --set-soname libgcc_s_f8.so libgcc_s_f8.so
# patchelf --set-soname libgfortran_f8.so libgfortran_f8.so
# patchelf --set-soname libquadmath_f8.so libquadmath_f8.so
# 
# # tell shared libs to look in '.' for other shared libs
# patchelf --set-rpath '$ORIGIN' libgcc_s_f8.so
# patchelf --set-rpath '$ORIGIN' libquadmath_f8.so
# patchelf --set-rpath '$ORIGIN' libgfortran_f8.so
# 
# # fix names to use local copies of shared libs
# patchelf --replace-needed libgcc_s.so.1 libgcc_s_f8.so libgfortran_f8.so
# patchelf --replace-needed libquadmath.so.0 libquadmath_f8.so libgfortran_f8.so
# 
# echo 'fixing Feff8 shared libraries'
# for F in libpotph.so libonepath.so libfeff8lpotph.so libfeff8lpath.so; do
#     patchelf --set-rpath '$ORIGIN' $F
#     patchelf --replace-needed libgcc_s.so.1 libgcc_s_f8.so $F
#     patchelf --replace-needed libquadmath.so.0 libquadmath_f8.so $F
#     patchelf --replace-needed libgfortran.so.3 libgfortran_f8.so $F
# done
# 
# echo 'fixing Feff8 binaries'
# for F in (feff6l feff8l_ff2x feff8l_genfmt feff8l_pathfinder
# 	  feff8l_pot feff8l_rdinp feff8l_xsph); do
#   patchelf --set-rpath '$ORIGIN' $F
#   patchelf --replace-needed libgcc_s.so.1 libgcc_s_f8.so $F
#   patchelf --replace-needed libquadmath.so.0 libquadmath_f8.so $F
#   patchelf --replace-needed libgfortran.so.3 libgfortran_f8.so $F
# 
# done
# 
