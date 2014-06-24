
from SCons.Environment import Environment
#from SCons.Script import *

def CompilationEnvironment():
    env = Environment()

    if env['FORTRAN'] == 'gfortran':
        env = Environment(FORTRANFLAGS = '-O3 -ffree-line-length-none')   ## -Wall -finit-local-zero
    elif env['FORTRAN'] == 'g77':
        env = Environment(FORTRANFLAGS = '-Wall -O2')
    elif env['FORTRAN'] == 'xlf':
        env = Environment(FORTRANFLAGS = '-qextern=trap')
    elif env['FORTRAN'] == 'ifort':
        env = Environment(FORTRANFLAGS = '-O3')

    return env
