## feff85exafs build system based on scons
## see HEADERS/license.h for license information

import SCons

from SCons.Environment import Environment
import sys
import os
if os.name == 'nt':
    from os import getcwd
else:
    from os import getcwd, chown
from os.path import realpath, join

from SCons.Script.SConscript import SConsEnvironment
if os.name != 'nt':
    from pwd import getpwnam, getpwuid
    from grp import getgrnam, getgrgid

DEBUG_ENV = False

## jsondir:
##   -I or -module flags: see
##   http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
##   n.b.: this gets evaluated in one of the subfolders, hence the ..

def CompilationEnvironment():
    """
    Determine how to build the Fortran parts of feff85exafs
    """
    args = {}
    jsondir = realpath(join(getcwd(), '..', 'json-fortran'))
    flags = {'jsondir': jsondir}

    gfortran_comp_flags = '-O3 -ffree-line-length-none -g -Wall'
    gfortran_link_flags = None

    json_comp_flags = ' -I{jsondir} -J{jsondir}'.format(**flags)
    json_link_flags = ' -L{jsondir} -ljsonfortran'.format(**flags)

    # darwin:
    if os.name == 'posix' and sys.platform=='darwin':
        args['platform'] = 'darwin'
        gfortran_comp_flags = '-O2 -arch x86_64 -Wall'
        gfortran_link_flags = '-dynamiclib -L/usr/local/gfortran/lib/ -lgfortran -lgfortranbegin'

    # windows: needs work!
    if os.name  == 'nt':
        #args['platform'] = 'Windows'
        pass

    env = Environment(**args)

    if env['FORTRAN'] == 'gfortran':
        # this was the suggestion in the top level Makefile in what
        # the FP gave us: "-O3 -ffree-line-length-none -finit-local-zero"
        # -O3 makes sense, as does -finit-local-zero.  I don't understand
        # the advantage of -ffree-line-length-none -- it seems to
        # encourage poor code style -- like we need more of that!
        # also,  -finit-local-zero fails on FMS/fmstot.f
        env = Environment(FORTRANFLAGS = gfortran_comp_flags + json_comp_flags, CFLAGS = '-g')
        if gfortran_link_flags is not None:
            env.Replace(SHLINKFLAGS = gfortran_link_flags + json_link_flags)
    elif env['FORTRAN'] == 'g77':
        env = Environment(FORTRANFLAGS = '-Wall -O2')
    elif env['FORTRAN'] == 'xlf':
        env = Environment(FORTRANFLAGS = '-qextern=trap')
    elif env['FORTRAN'] == 'ifort':
        ## I think the -module flg is correct ... untested ...
        env = Environment(FORTRANFLAGS = '-O3 -module '+jsondir)

    if DEBUG_ENV:
        for key, val in env.items():
            try:
                print( key, val)
            except:
                pass
        sys.exit()

    if os.name == 'nt':
        env.PrependENVPath('PATH', os.environ['PATH'])
    
    return env

## need to be able to get prefix from command line
def InstallEnvironment():
    """
    Determine installation locations for libraries and executables.
    """
    ienv = Environment()
    #prefix = ARGUMENTS.get('prefix', '/usr/local')
    if os.name == 'nt':
        import larch
        prefix = larch.larchlib.sys_larchdir
        dlldir = larch.larchlib.get_dlldir()
        # Here are our installation paths:
        ienv['i_prefix'] = prefix
        ienv['i_lib']    = join(prefix, 'dlls', dlldir)
        ienv['i_bin']    = join(prefix, 'bin')
        ienv['i_inc']    = join(prefix, 'include')
        ienv['i_data']   = join(prefix, 'share')
    else:
        prefix = '/usr/local'
        # Here are our installation paths:
        ienv['i_prefix'] = prefix
        ienv['i_lib']    = join(prefix, 'lib')
        ienv['i_bin']    = join(prefix, 'bin')
        ienv['i_inc']    = join(prefix, 'include')
        ienv['i_data']   = join(prefix, 'share')
    return ienv


def FindOtherObjects(deplist, env, kind):
    objsuff = env['SHOBJSUFFIX']
    if kind == 'static': objsuff = env['OBJSUFFIX']
    out = []
    for dep in deplist:
        dname, fname = dep.split('/')
        out.append(join('..', dname, fname + objsuff))
    return out

##+----------------------------------------------------------------------------------------------------
## This implementation of a Chown factory closely follows Chmod from /usr/lib/scons/SCons/Defaults.py
## and the discussion from http://www.scons.org/wiki/InstallTargets

def get_paths_str(dest):
    # If dest is a list, we need to manually call str() on each element
    if SCons.Util.is_List(dest):
        elem_strs = []
        for element in dest:
            elem_strs.append('"' + str(element) + '"')
        return '[' + ', '.join(elem_strs) + ']'
    else:
        return '"' + str(dest) + '"'

## FIXME.MAYBE: recognize numeric and string uid and gid
def chown_func(dest, owner, group):
    SCons.Node.FS.invalidate_node_memos(dest)
    if not SCons.Util.is_List(dest):
        dest = [dest]
    for element in dest:
        #os.chmod(str(element), getpwnam(owner).pw_uid, getgrnam(owner).gr_gid)
        chown(str(element), owner, owner)

def chown_strfunc(dest, owner, group):
    return 'Chown(%s, %s.%s)' % (get_paths_str(dest), getpwuid(owner).pw_name, getgrgid(group).gr_name)

SConsEnvironment.Chown = SCons.Action.ActionFactory(chown_func, chown_strfunc)

def InstallOwner(env, dest, files, owner, group):
    """
    Used for files to be owned by the normal user even if the installation is run as sudo
    """
    obj = env.Install(dest, files)
    for i in obj:
        env.AddPostAction(i, env.Chown(str(i), owner, group))
    return dest

SConsEnvironment.InstallOwner = InstallOwner
