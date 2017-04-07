Materials for unit testing
==========================

For each example material, several files are provided:

 1. Structural data, either in the form of an
    [Atoms input file](http://bruceravel.github.io/demeter/artug/feff/index.html#crystaldata),
    a [CIF file](http://www.iucr.org/resources/cif), or a
    [Feff input file](http://monalisa.phys.washington.edu/feff/wiki/static/f/e/f/FEFFinp_4993.html).

 2. Where crystal data is given in the form of an Atoms input or CIF
	file, the crystal data has been converted into a file very much
	like a Feff input file and with the extension `.mustache`.

 3. The Feff input files are not yet in a form ready for running Feff.
    Each one is in a form used by the
    [Mustache](http://mustache.github.io/) templating system.  A few
    bits of input data have been replaced by tokens that look like
    this: `{{scf}}`.  These tokens will be replaced by values
    appropriate to the test being performed.  Those values are usually
    taken from a `.json` file with the same name as the folder itself.
    The mustache files here were edited by hand or made using
    [the `feff85exafs` template in Demeter](https://github.com/bruceravel/demeter/blob/master/lib/Demeter/templates/atoms/feff85test.tmpl).

 4. For tests which have data associated with them and will include
    results of fits to EXAFS data as part of their test, the data is
    provided in the form of an
    [Athena project file](http://bruceravel.github.io/demeter/aug/output/project.html)
    and as a column ASCII data containing chi(k).  The column data
    file is in the
    [XDI format](https://github.com/XraySpectroscopy/XAS-Data-Interchange).

 5. For tests which have data associated with them and will include
    results of fits to EXAFS data as part of their test, a file of
    python code defining the fit must be provided.  This **must**
    define a function called `do_fit` which implements the fit in
    larch and provides a plot and fit report for interactive use.

 6. Most of the material folders also contain some kind of image
    showing the nature of the coordination environment about the
    absorbing atom.

 7. Each materials folder has a subfolder containing the baseline Feff
    calculation.  This is a run of Feff 8.5 as delivered by the Feff
    Project.  The sense in which this is a baseline is that, as
    changes are made to the Feff85EXAFS code base, those changes can
    be tested against this original state of the software.  Changes to
    the code base should not result in changes to the output of Feff.
    This baseline can also serve as a platform for testing changes
    across versions of Feff, allowing us to probe in a systematic way
    the differences in EXAFS analysis between versions 6 through 9 or
    Feff.  Baseline calculations have been made with and without
    self-consistency.


The testing infrastructure is implemented with
[Nose](https://nose.readthedocs.org/en/latest/index.html) and a
[Larch](https://github.com/xraypy/xraylarch) plugin.

The testing infrastructure is designed so that it is easy to add new
tests, so long as an appropriate set of files is provided for each new
test.

Note that the file naming conventions for the test files are quite
strict.  If you introduce a new material, say Ceria (CeO2), and you
name the folder containing its files `Ceria/`, then the you must
follow these rules:

1. The structure file **must** be `Ceria.<extension>` or
   `Ceria_atoms.inp`.  Here `<extension>` is something like `cif` or
   `sdf`.
2. The Mustache template file **must** be called `Ceria.mustache`.
3. The chi(k) **must** be called `Ceria.chik`, the first data column
   **must** be wavenumber and the second column **must** be
   un-k-weighted chi(k).  You **should** save the chi(k) file as an
   XDI file. 
4. The JSON file with the configuration for the Feff run **must** be
   called `Ceria.json`.
5. There **must** be a folder called `baseline/` which has folders
   called `withSCF/` and `noSCF/`. These contain the baseline
   calculations with and without self-consistency.  The baseline
   calculation should be made using a point in the feff85exafs history
   ([this point, for example](https://github.com/xraypy/feff85exafs/commit/cac0f8c90749ce52581a658c5a6c8ae144cc2211))
   from before changes were made to the code as it was delivered to us
   from the Feff Project.
6. If you provide an Athena project file, it should be called
   `Ceria.prj`.
7. If fits to data are part of the test, there must be a file called
   `Ceria.py` containing python code defining the fit to the data.
8. You must add `Ceria` to the `folders` tuple at the top of
   `tests/test_materials.py`.
9. You should provide a `README.md` file with basic information in
   markdown format.  Any other files, for instance images displayed in
   the `README.md` file can have any name (since they will not be used
   in testing)

---

## Materials:

1. **copper metal**: because copper

2. **nickel oxide, NiO**: this is a simple, cubic, metal oxide.  It
   represents a problem slightly more complicated than copper metal.
   It also is easy to get a decent fit out to the sixth coordination
   shell.

3. **uraninite, UO2**: this is an f-electron system

4. **zircon, ZrSiO4**: this has Si, a tender energy absorber, with a
   4d backscatterer

5. **bromoadamantane**: this is a small molecule which can be fit with a
   fairly simple model of four paths, but for which the 6 nearby hydrogen
   scatterers play a big role in the fit.

6. **ferrocene macrocycle**: this is a iron organometallic with Fe in
   between two 5 member carbon rings, so it has 10 C neighbors at a range
   of distances from 2.026 A to 2.099 A

7. **lanthanum cuprate**: this is an orthorhombic crystal with the
   copper atom in an octahedral configuration.  The oxygen atoms along
   the z axis are much farther away than the ones in the plane.  This,
   then, is a test that involves Feff's polarization and ellipticity
   calculations.  There is one folder for calculating with the
   incident light parallel to the z axis (polarization only) and one
   for light perpendicular to the z axis (polarization and
   ellipticity).

---

## Installing and using the unit testing tool

Make sure that all parts of Feff have been compiled successfully.  The
unit test framework currently uses larch's `feffrunner` class to make
the test runs of Feff.  You must have the termcolor, pystache and nose
python libraries installed.  (In debian/ubuntu, these are called
`python-termcolor`, `python-pystache`, and `python-nose`.  Or you can
use [pip](https://pip.pypa.io/en/stable/)).

(Eventually, this testing framework will need to test the data tables
from the wrapper against the `feffNNNN.dat` files from the baseline
calculation.)

Copy the file `f85ut.py` to the larch plugins folder
(`$HOME/.larch/plugins/` on Unix; `C:\Users\<ME>\larch\plugins` or
`C:\Program Files\larch\plugins` on Windows).  Alternately, make a
symbolic link from the plugin folder to the file:

      ln -s  `pwd`/f85ut.py ~/.larch/plugins/f85ut.py

Here, I have assumed that you will run this command in the folder
containing both this README file and the `f85ut.py` file.  That is why
the call to `pwd` works.



---

### Run tests through nose

In the `materials` folder, run `nosetests` at the command line.  Nose
writes a report on the results of the test sequence.

Any tests that fail can be further examined interactively within
Larch.

When run through Nose, the beginning of the test sequence is *very*
time consuming as all the Feff calculations are made before any of the
actual tests are made.  I find it helpful to run `nosetests
--verbosity=3`, which gives some feedback about what is actually
happening.


#### Run Feff with self-consistency

set the `FEFF_TEST_SCF` environment variable to `True`.  With bash,
zsh, etc:

    ~> export FEFF_TEST_SCF=True

With csh, tcsh, etc:

    ~> setenv FEFF_TEST_SCF "True"

#### Skip canned data tests

Normally, a data test will be run if the `<name>.py` file is present
in the material's folder.  The data test for that material will be
skipped if a file called `<name>.skip` exists in its folder.  For
example, to skip the test for UO2, do

	~> touch UO2/UO2.skip

To re-enable the data test for UO2, do

	~> rm UO2/UO2.skip



---

### Interactive testing in Larch

Fire up Larch from the folder containing this README file and do the
following:

    larch> add_plugin('f85ut')
    larch> my_ut = ut('Copper')

That will create a unit testing object called `my_ut` and point it at
the folder containing the materials related to unit tests on a
calculation for copper metal.

Next do

     larch> my_ut.run()

This will run a calculation on copper metal using Feff85exafs in its
current state.  Once this is done, you can check the results of the
calculation against the baseline calculation, i.e. a calculation made
using Feff85exafs as it was delivered to us by the Feff Project.

To run Feff with self-consistency, do

	 larch> my_ut.doscf=True

The Feff run is, of course, much more time consuming with
self-consistency.

To see whether the calculation of a path differs from the baseline
calculation, do

     larch> my_ut.compare(path_index)

where `path_index` is the `NNNN` in `feffNNNN.dat`.

This will compute a simple R factor between the magnitude of F\_eff in
the baseline and the test run.  It will also compute and R factor for
the phase of F\_eff.  If `my_ut.doplot` is True (which is the default),
plots of magnitude and phase of F\_eff will be made including both the
baseline and the test run.

You can compare other parts of the calculation:

     larch> my_ut.compare(path_index, part)

where `part` is one of

part    | purpose
--------|------------------------------------------------
feff    | test the magnitude AND phase of F_eff (default)
amp     | test the total amplitude of the path
phase   | test the total phase of the path
lambda  | test the mean free path
caps    | test the central atom phase shift
redfact | test the reduction factor
rep     | test the real part of the complex wavenumber

To test if a path index was saved from the Feff calculation

     larch> if my_ut.available(nnnn):
     larch>     my_ut.compare(nnnn)
     larch> end if

To examine various other quantities from the Feff calculation:

     larch> my_ut.feffterms()
	 larch> print my_ut.radii('testrun', 'muffintin') my_ut.radii('baseline', 'muffintin') 
	 larch> print my_ut.radii('testrun', 'norman')    my_ut.radii('baseline', 'norman') 
	 larch> print my_ut.s02('testrun')                my_ut.s02('baseline') 

Some of the materials have data tests.  This

     larch> my_ut.fit()

runs a canned fit, once using the baseline Feff calculation and once
using the test run.  You can then compare fitting parameters and
statistics from the two.  The fit groups are `my_ut.blfit` and
`my_ut.trfit`.  You could, for example, examine the `amp` parameter:

     larch> print my_ut.blfit.params.amp.value my_ut.blfit.params.amp.stderr
     larch> print my_ut.trfit.params.amp.value my_ut.trfit.params.amp.stderr

Finally, clean up the test run by doing:

     larch> my_ut.clean()

Some convenience functions exported by the plugin:

* use `my_ut = ir('Copper')` to define the unit testing object and run feff
* use `my_ut = irc('Copper')` to define the unit testing object, run
  feff, and make a comparison for first path
* use `my_ut = irf('Copper')` to define the unit testing object, run
  feff, and run the fit


# Still to do

* use wrapper to do tests of baseline feffNNNN.dat file to in-memory
  data table

* capture and interpret Feff's screen messages to use number of SCF
  iterations as a unit test

* More materials:
	+ A polymer, i.e. something pseudo-one-dimensional
	+ Something from the first row of the periodic table
	+ Something with lead as the absorber
	+ Something from column 1 or column 2 of the periodic table
	+ Something with a ring structure and strong, high-order MS paths
      (paradibromobenzene, perhaps)
	+ americium, a transuranic that was above Feff6's Z cutoff
