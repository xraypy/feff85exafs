List of materials for unit testing
==================================

For each example material, several files are provided:

 1. structural data, either in the form of an
    [Atoms input file](http://bruceravel.github.io/demeter/artug/feff/index.html#crystaldata),
    a [CIF file](http://www.iucr.org/resources/cif), or a
    [Feff input file](http://monalisa.phys.washington.edu/feff/wiki/static/f/e/f/FEFFinp_4993.html).

 2. Where crystal data is given in the form of an Atoms input or CIF
    file, the crystal data has been converted into a Feff input file.

 3. The Feff input files are not yet in a form ready for running feff.
    Each one is in a form used by the
    [Mustache](http://mustache.github.io/) templating system.  A few
    bits of input data have been replaced by tokens that look like
    this: `{{scf}}`.  These tokens will be replaced by values
    appropriate to the test being performed.  Those values are taken
    from a `.json` file with the same name as the folder itself.

 4. For test which have data associated with them and, so, will
    included results of fits to EXAFS data as part of their test, the
    data is provided in the form of an
    [Athena project file](http://bruceravel.github.io/demeter/aug/output/project.html)
    and as a column ASCII data containing chi(k).  The column data
    file is in the
    [XDI format](https://github.com/XraySpectroscopy/XAS-Data-Interchange).

 5. Most of the material folders also contain some kind of image
    showing the nature of the coordination environment about the absorbing atom.

 6. Each materials folder has a folder containing the baseline Feff
    calculation.  This is a run of Feff 8.5 as delivered by the Feff
    Project.  The sense in which this is a baseline is that, as
    changes are made to the Feff85EXAFS code base, those changes can
    be tested against this original state of the software.  Changes to
    the code base should not result in changes to output of Feff.
    This baseline can also serve as a platform for testing changes
    across versions of Feff, allowing us to probe in a systematic way
    the differences in EXAFS analysis between versions 6 through 9 or
    Feff.


The testing infratructure is implemented as a
[Larch](https://github.com/xraypy/xraylarch) plugin.  See below for
details.

The testing infratructure will be designed so that it is easy to add
new tests, so long as an appropriate set of files is provided for each
new test.

Please note that the file naming comventions for the test files are
quite strict.  If you introduce a new material, say Ceria (CeO2), and
you name the folder containing its files either `Ceria`, then the
following must be true:

1. The structure file **must** be `Ceria.<extension>` or
   `Ceria_atoms.inp`.  Here `<extension>` is something like `cif` or
   `sdf`.
2. The Mustache template file **must** be called `Ceria.mustache`.
3. The chi(k) **must** be called `Ceria.chik`, the first data column
   **must** be wavenumber and the second column **must** be
   un-k-weighted chi(k).
4. The JSON file with the configuration for the Feff run **must** be
   called `Ceria.json`.
5. If you provide an Athena project file, is should be called
   Ceria.prj.
6. You should provide a `README.md` file with basic information in
   markdown format.  Any other files, for instance images displayed in
   the `README.md` file can have any name (since they will not be used
   in testing)

## Materials:

1. **copper metal**: because copper

2. **nickel oxide, NiO**: this is a simple, cubic, metal oxide.  It
   represents a problem slightly more complicated than copper metal.

3. **uraninite, UO2**: this is an f-electron system

4. **zircon, ZrSiO4**: this has Si, a tender energy absorber, with a 4d backscatterer

5. **bromoadamantane**: this is a small molecule which can be fit with a
   fairly simple model of four paths, but for which the 6 nearby hydrogen
   scatterers play a big role in the fit.

6. **ferrocene macrocycle**: this is a iron organometallic with Fe in
   between two 5 member carbon rings, so it has 10 C neighbors at a range
   of distances from 2.026 A to 2.099 A

## Installing and using the unit testing tool

Copy the file `f85ut.py` to the larch plugins folder (either
`$HOME/.larch/plugins/` or `/usr/local/share/larch/plugins` on Unix,
or `C:\Users\ME\larch\plugins` or `C:\Program Files\larch\plugins` on
Windows).  Alternately, make a symbolic link from the plugin folder to
the file:

      ln -s  `pwd`/f85ut.py ~/.larch/plugins/f85ut.py

Here, I have assumed that you have run this command in the folder
containing both this README file and the `f85ut.py` file.  That is
what lets the call to `pwd` work.

Fire up Larch and do the following:

    add_plugin('f85ut')
    my_ut = ut('Copper')

That will creat a unit testing object called `my_ut` and point it at
the folder containing the materials related to unit tests on a
calculation for copper metal.

Next do

     my_ut.run()

This will run a calculation on copper metal using Feff85exafs in its
current state.  Once this is done, you can check the results of the
calculation against the baseline calculation, i.e. a calculation made
using Feff85exafs as it was delivered to us by the Feff Project.

To see whether the calculation of the first path differs from teh
baseline calculation, do

     my_ut.compare(1)

This will compute a simple R factor between the magnitude of F\_eff in
the baseline and the test run.  It will also compute and R factor for
the phase of F\_eff.  If my_ut.doplot is True (which is the default),
plots of magnitude and phase of F\_eff will be made including both the
baseline and the test run.

You can compare other parts of the calculation:

     my_ut.compare(1, part)

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

You can test to see if a path index was saved from the Feff calculation

      if my_ut.available(nnnn):
	      my_ut.compare(nnnn)

Finally, clean up the test run by doing:

      my_ut.clean()

