
from SCons.Environment import Environment

def CompilationEnvironment():
    env = Environment()

    if env['FORTRAN'] == 'gfortran':
        # this was the suggestion in the top level Makefile in what
        # the FP gave us: "-O3 -ffree-line-length-none -finit-local-zero"
        # -O3 makes sense, as does -finit-local-zero.  I don't understand
        # the advantage of -ffree-line-length-none -- it seems to 
        # encourage poor code style -- like we need more of that!
        # also,  -finit-local-zero fails on FMS/fmstot.f
        env = Environment(FORTRANFLAGS = '-O3 -ffree-line-length-none -Wall')   ## -Wall -finit-local-zero
    elif env['FORTRAN'] == 'g77':
        env = Environment(FORTRANFLAGS = '-Wall -O2')
    elif env['FORTRAN'] == 'xlf':
        env = Environment(FORTRANFLAGS = '-qextern=trap')
    elif env['FORTRAN'] == 'ifort':
        env = Environment(FORTRANFLAGS = '-O3')

    return env
