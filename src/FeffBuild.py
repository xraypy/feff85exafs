## feff85exafs build system based on scons
## see HEADERS/license.h for license information

from SCons.Environment import Environment
from os      import getcwd
from os.path import realpath, join

## jsondir:
##   -I or -module flags: see 
##   http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
##   n.b.: this gets evaluated in one of the subfolders, hence the ..

def CompilationEnvironment():
    """
    Determine how to build the Fortran parts of feff85exafs
    """
    env = Environment()
    jsondir = realpath(join(getcwd(), '..', 'JSON'))

    if env['FORTRAN'] == 'gfortran':
        # this was the suggestion in the top level Makefile in what
        # the FP gave us: "-O3 -ffree-line-length-none -finit-local-zero"
        # -O3 makes sense, as does -finit-local-zero.  I don't understand
        # the advantage of -ffree-line-length-none -- it seems to 
        # encourage poor code style -- like we need more of that!
        # also,  -finit-local-zero fails on FMS/fmstot.f
        env = Environment(FORTRANFLAGS = '-O3 -ffree-line-length-none -Wall -I'+jsondir)   ## -pedantic -finit-local-zero
    elif env['FORTRAN'] == 'g77':
        env = Environment(FORTRANFLAGS = '-Wall -O2')
    elif env['FORTRAN'] == 'xlf':
        env = Environment(FORTRANFLAGS = '-qextern=trap')
    elif env['FORTRAN'] == 'ifort':
        ## I think the -module flg is correct ... untested ...
        env = Environment(FORTRANFLAGS = '-O3 -module '+jsondir)

    return env

## need to be able to get prefix from command line
def InstallEnvironment():
    """
    Determine installation locations for libraries and executables.
    """
    ienv = Environment()
    #prefix = ARGUMENTS.get('prefix', '/usr/local')
    prefix = '/usr/local'
    # Here are our installation paths:
    ienv['i_prefix'] = prefix
    ienv['i_lib']    = prefix + '/lib'
    ienv['i_bin']    = prefix + '/bin'
    ienv['i_inc']    = prefix + '/include'
    ienv['i_data']   = prefix + '/share'
    return ienv
