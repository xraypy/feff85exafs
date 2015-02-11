## source prep.sh

alias ls='ls -F --color=always'
alias scons='python /c/Python27/Scripts/scons.py'
alias rm='rm -i'

## need to put the location of gcc 4.9 at the front of the path
## (Demeter, for instance, comes with 4.7, which won't work for json-fortran)
export PATH=/C/Program\ Files/mingw-w64/x86_64-4.9.2-win32-seh-rt_v3-rev1/mingw64/bin:$PATH

## want the installation location of swig at the end of the path
export PATH=$PATH:/c/Program\ Files/swigwin-3.0.5
export PATH=$PATH:/c/Program\ Files/larch/bin
export PATH=$PATH:/c/Program\ Files/larch/dlls/win32

export VS90COMNTOOLS="/C/Users/bravel/AppData/Local/Programs/Common/Microsoft/Visual C++ for Python"

## these need to be in the path as well
# /c/Python27
# /c/Program Files/Larch/bin
# /c/Program Files/swigwin-3.0.5
