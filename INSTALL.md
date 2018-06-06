# Feff8L installation

Building and Installing  Feff8L, aka feff85exafs:


## Requirements

* gfortran v 4.9 or higher (someone please verify) and an accompanying gcc compiler.

* Python, version 2.7 or 3.5 or higher.  Some scripts may still be
  incompatible with Python3 -- please let us know if you run into a problem
  with this.

## Download source code

Download the tarball (or zip ball) for this project, and unpack onto your
computer,  or use `git clone` on the master branch of this repository to
pull onto your computer. You may be able to use one of the following
commands (you do not need all of them!):

usig wget:

```
    wget  https://github.com/xraypy/feff85exafs/archive/0.2.tar.gz
    tar xvf 0.2.tar.gz
    cd feff85exafs-0.2
```

using curl:

```
    curl -L https://github.com/xraypy/feff85exafs/archive/0.2.tar.gz
    tar xvf 0.2.tar.gz
    cd feff85exafs-0.2
```

using git:

```
     git clone https://github.com/xraypy/feff85exafs.git
     cd feff85exafs
```

## Build and install

It should be possible to build feff8l with

```
    make install
```

By default, this will place the executables in the `local_install` folder
in the `feff85exafs` folder.  You should be able to install these elsewhere
simply by copying them to your favorite location of binaries, perhaps with

```
    sudo cp -pr local_install/bin/* /usr/local/bin/.
    sudo cp -pr local_install/lib/* /usr/local/lib/.
    sudo cp -pr local_install/include/* /usr/local/include/.
```

You can also change the install location by editing the first few lines of the
Makefile (in this directory).  By default this reads

```
    export PREFIX  = /usr/local
    export PREFIX =  ${CURDIR}/local_install
```

But you can change this to just read

```
    export PREFIX  = /usr/local
```

or anywhere else you would like to install to.  The executables will go in
the `bin/` folder below this main install location.

For further details on compiling Feff8l and using
[json-fortran](https://github.com/jacobwilliams/json-fortran),
see [`src/README.md`](src/README.md).


## Testing

Testing requires
[python's nose tool](https://nose.readthedocs.org/en/latest/) and [Larch](https://github.com/xraypy/xraylarch).

To compile and  test, do the following in the top directory

```
	  feff8l> make
	  feff8l> LD_LIBRARY_PATH='wrappers/python:src/GENFMT' && nosetests --verbosity=3
	  feff8l> sudo make install
```


For more example programs using the fortran entry point to the stand-alone
F_eff calculations or programs using the C wrapper see the `wrappers/`
directory.  There you will also find language bindings to the C wrapper,
including python and perl.

# Running Feff8L


With feff8l so that it is in your path, you can create a folder and place a
`feff.inp` folder in that folder, and then run the `feff8l` python script
on that directory:

```
      mkdir FeffTest
      cp SomeFeff8.inp FeffTest/.
      feff8l  FeffTest/.
```

feff8l is a Python 2 script that runs all (or selected) modules of Feff85L in order:
     feff8l_rdinp
     feff8l_pot
     feff8l_xsph
     feff8l_pathfinder
     feff8l_genfmt
     feff8l_ff2x

After running the feff8l script, the directory FeffTest should contain the
Feff8L output files. Since they have generic names, it is recommended that
you copy them into a sensibly named subdirectory in your EXAFS Project
Directoy.

To understand the feff.inp file, have a look at feff6l_doc.txt or feff8.tex in the directory feff85exafs/doc/ .
