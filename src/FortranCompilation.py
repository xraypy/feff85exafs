## feff85exafs build system based on scons
## see HEADERS/license.h for license information

from SCons.Environment import Environment
from os      import getcwd
from os.path import realpath, join

## -I or -module flags: see 
## http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
jsondir = realpath(join(getcwd(), '..', 'JSON'))
## n.b.: this gets evaluated in one of the subfolders, hence the ..

#prefix   = ARGUMENTS.get('prefix', '/usr/local')

def CompilationEnvironment():
    env = Environment()

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

    # Here are our installation paths:
    # env['i_prefix'] = prefix
    # env['i_lib']    = prefix + '/lib'
    # env['i_bin']    = prefix + '/bin'
    # env['i_inc']    = prefix + '/include'
    # env['i_data']   = prefix + '/share'

    return env
