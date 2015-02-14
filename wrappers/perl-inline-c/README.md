# NAME

Xray::FeffPath - Moose wrapper around an Inline::C interface to libfeffpath

# VERSION

1.00

# SYNOPSIS

This provides a simple object oriented interface to the `feffpath`
wrapper around Feff's `onepath` subroutine.  `onepath` combines the
F-matrix calculation in `GENFMT/genfmt.f` with the construction of the
common presentation of F-effective from `FF2X/feffdt.f`.

Given the geometry of a scattering path in the form of Cartesian
coordinates of atoms making up the path, along with some additional
information such as the degeneracy, compute F-effective for that
scattering path.

This replaces the need to read a `feffNNNN.dat` file.  Given a
scattering geometry, the equivalent data can be computed directly.

The following computes the columns of the `feff0004.dat` file for
copper metal:

    #!/usr/bin/perl
    use strict;
    use warnings;
    use Xray::FeffPath;

    my $path = Xray::FeffPath->new();
    $path->degen(48);
    $path->Index(4);
    $path->phpad('../fortran/phase.pad');
    $path->atom(0, 0, -3.61, 1);
    $path->atom(-1.805, 0, -1.805, 1);
    $path->path;
    my @kgrid  = @{ $path->k };
    my @caps   = @{ $path->real_phc };
    my @amff   = @{ $path->mag_feff };
    my @phff   = @{ $path->pha_feff };
    my @redfac = @{ $path->red_fact };
    my @lambda = @{ $path->lam };
    my @rep    = @{ $path->rep };
    undef $path

# INSTALLATION

After you have built and installed _feff85exafs_, do the following:

    perl Makefile.PL
    make
    make test
    sudo make install

That's it!  Note, though, that building this wrapper **requires** that
the fortran and C code _feff85exafs_ be completely compiled and that
the resulting libraries (and other files) be successfully installed.

# METHODS

- `new`

    Create the FeffPath object

        my $path = Xray::FeffPath->new();

- `atom`

    Add an atom to the scattering path, supplying the Cartesian
    coordiantes of the atom and its unique potential index.

        $nleg = $path->atom($x, $y, $z, $ipot);

    This returns the current number of legs in the path.  To make a SS
    paths do:

        $path->atom($x,  $y,  $z, $ipot1);
        print $path->nleg, $/;
        #   ==prints===> 2

    To make a MS path do:

        $path->atom($x1, $y1, $z1, $ipot1);
        $path->atom($x2, $y2, $z2, $ipot2);
        $path->atom($x3, $y3, $z3, $ipot3);
        print $path->nleg, $/;
        #   ==prints===> 4

    Note that there is no setter method for the nleg attribute.  Rather,
    `nleg` is incremented by each call to the `atom` method.

- `path`

    Compute the scattering path

        $path->path;

- `clear`

    Reinitialize the scattering path

        $path->clear;

To destroy a FeffPath object, simply

    undef $path;

# ATTRIBUTES

Each of the members of the FEFFPATH struct is wrapped in a getter
method of the same name.  Each of the following exists:

## Scalar valued attributes

- `phpad` (character, default = phase.pad)

    The path to the `phase.pad` file.  This gets right-padded with spaces
    to 256 characters, which is the length of the parameter in the
    underlying Fortran library.  Therefore, the value returned by the
    getter will be right-padded.

- `Index` (integer, default = 9999)

    The path index, used for form the name of the output `feffNNNN.dat`
    file.

    Accessing this and other scalar-values attributes works in the typical
    Moose-y way:

        $path->Index(4);
        print $path->Index, $/;
        ## ==prints==> 4

    Note that this attribute (but only this one) is capitalized to avoid
    confusion with the built-in `index` function.

- `degen` (float, default = 0)

    The path degeneracy.  This is required input for a calculation.

- `nleg` (integer, default = 0)

    The number of legs in the path.  This should never be set by hand.  It
    gets incremented by calls to the `atom` method.

- `iorder` (integer, default = 2)

    Order of approximation for multiple scattering paths.

- `nnnn` (boolean, default = 0)

    Flag for writing `feffNNNN.dat` file.

- `json` (boolean, default = 0)

    Flag for writing `feffNNNN.json` file.

- `verbose` (integer, default = 0)

    Flag for writing screen messages as files are written.

- `ipol` (boolean, default = 0)

    Perform polarization calculation.  This gets toggled to true whenever
    any of `elpty`, `evec`, or `xivec` are set.

- `elpty` (integer, default = 0)

    Ellipticity for use in polarization calculation.

- `reff` (float, default = 0)

    Half path length.  This should never be set by hand.  It
    gets set correctly when the `path` method is called.

- `ne` (integer, default = 0)

    Number of outpit k grid point.  This should never be set by hand.  It
    gets set correctly when the `path` method is called.

## Error reporting

- `errorcode` (integer)

    An integer error code from either `atom` or `path`.  Check for
    non-zero value.  For a non-zero return from `atom`, calling `path`
    is likely to cause a problem in the genfmt calculation.  For a
    non-zero return from `path`, the genfmt calculation was skipped.  See
    `errormessage` for an explanation of the problem.

- `errormessage` (string)

    A explanation of words of the problem found during `atom` or `path`.

        $path->degen(-48);
        $path->atom(0, 0, 3.61, 1);
        $path->path;
        if ($path->errorcode) {
           print $path->errormessage;
        };
        ## ==prints==>
         Error in make_path
             path degeneracy (-48.00) is negative

## Attributes related to the potential model

Several parameters are captured by the C struct (and therefore
available via the perl wrapper) which are related to how Feff's
potential model was calculated.  This information is written to the
`feffNNNN.dat` header.

- `version`

    A string identifying the version number and the \`feffpath\` wrapper
    revision.

- `exch`

    A string identifying the type of potential model used.  The
    possibilities are

        H-L exch (the default)
        D-H exch
        Gd state
        DH - HL
        DH + HL
        val=s+d
        sigmd(r)
        sigmd=c

    See the description of the EXCHANGE keyword in the Feff document
    [http://feffproject.org/feff/Docs/feff9/feff90/feff90\_users\_guide.pdf](http://feffproject.org/feff/Docs/feff9/feff90/feff90_users_guide.pdf)
    at page 92.

- `edge`

    A poor estimate of the energy threshold relative to the atomic value, in eV.

- `gam_ch`

    The core-hole lifetime, in eV.

- `kf`

    The k value at the Fermi energy, in inverse Angstroms.

- `mu`

    The Fermi energy, in eV.

- `rnorman`

    The average Norman radius.

- `rs_int`

    An estimate of the interstitial radius.

- `vint`

    The interstitial potential, in eV.

## Constants from feffpath.h

The following constants are captured from the feffpath header file:

- `nex` (150)

    The maximum number of points in the energy grid.

- `nphx` (11)

    The maximum number of unique potentials.

- `npatx` (8)

    The maximum number of path atoms

- `legtot` (9)

    The maximum number of legs in a path (npatx + 1)

- `bohr` (0.529177249)

    A unit of length, in Angstrom.

- `ryd` (13.605698)

    A Rydberg, a unit of energy, in eV.

- `hart` (27.211396)

    A Hartree, a unit of energy, in eV.

## Attributes with values of array-reference

- `absorber`

    The Cartesian coordinates of the absober atom.  By default the
    absorber is placed at (0,0,0).

- `ipot`

    A list of the unique potentials of the scatterers in a path.  This is
    populated by the calls to the `atom` method.

- `rat`

    A list of lists of the Cartesian coordinates of the scatterers in a
    path.  This is populated by the calls to the `atom` method.

- `evec`

    The electric vector of the incident beam for the polarization
    calculation.  Note that the setter's argument is an array reference as
    is the getter's return value.

        $path->evec([0,0,1]);
        $evec_ref = $path->evec;
        print join(", ", @$evec_ref);
        ## ==prints==> 0, 0, 1

    Setting `evec` will also set `ipol` to true.

- `xivec`

    The Poynting vector of the incident beam for the ellipticity
    calculation.  Note that the setter's argument is an array reference as
    is the getter's return value.

        $path->xivec([1,1,0]);
        $xivec_ref = $path->xivec;
        print join(", ", @$xivec_ref);
        ## ==prints==> 1, 1, 0

    Setting `xivec` will also set `ipol` to true.

- `ri`

    The leg lengths (in Angstroms) of the specified path.  This gets set
    when `path` is called.  It returns an reference to an array with
    `nleg` elements.

- `beta`

    The Eulerian beta angles (in degrees) of the specified path.  This
    gets set when `path` is called.  It returns a reference to an array
    with `nleg` elements.

- `eta`

    The Eulerian eta angles (in degrees) of the specified path.  This gets
    set when `path` is called.  It returns a reference to an array with
    `nleg` elements.

- `k`

    A reference to an array containing the k grid of the calculation.
    This is column 1 from `feffNNNN.dat`.

- `real_phc`

    A reference to an array containing the central atom phase shifts.
    This is column 2 from `feffNNNN.dat`.

- `mag_feff`

    A reference to an array containing the magnitude of F-effective.  This
    is column 3 from `feffNNNN.dat`.

- `pha_feff`

    A reference to an array containing the phase of F-effective.  This is
    column 4 from `feffNNNN.dat`.

- `red_fact`

    A reference to an array containing the reduction factor.  This is
    column 5 from `feffNNNN.dat`.

- `lam`

    A reference to an array containing lambda, the mean free path.  This is
    column 6 from `feffNNNN.dat`.

- `rep`

    A reference to an array containing the real part of the complex
    momentum.  This is column 7 from `feffNNNN.dat`.

# Error codes

The error codes returned by the `atom` and `path` methods are the
sums of the codes actually found by the method call, thus the returned
error codes are meant to be interpreted bitwise.

## `atom` method

If any of these are triggered, it is very likely that a call to
`path` will return unreliable results or may crash the program.

- _1_

    ipot argument to add\_scatterer is less than 0

- _2_

    ipot argument to add\_scatterer is greater than 7

- _4_

    coordinates are for an atom too close to the previous atom in the path

- _8_

    nlegs greater than legtot

## `path` method

If any of these are triggered, genfmt will not be called and all
output arrays will be zero-filled.

- _1_

    the first atom specified is the absorber

- _2_

    the last atom specified is the absorber

- _4_

    path degeneracy is negative

- _8_

    path index not between 0 and 9999

- _16_

    ellipticity not between 0 and 1

- _32_

    iorder not between 0 and 10

- _64_

    `phase.pad` cannot be found or cannot be read

# EXTERNAL DEPENDENCIES

- [Inline](https://metacpan.org/pod/Inline), and [Inline::C](https://metacpan.org/pod/Inline::C)
- [Moose](https://metacpan.org/pod/Moose)
- [MooseX::NonMoose](https://metacpan.org/pod/MooseX::NonMoose) and [MooseX::Aliases](https://metacpan.org/pod/MooseX::Aliases)

# BUGS AND LIMITATIONS

- Setting a boolean unsets phpad at Inline::C level.  Weird!
- Polarization test is quite trivial.  Should test against an actual
calculation.
- `ipol` and `rat` not accessed via wrapper.  That said,
Xray::FeffPath keeps those attributes current.

Please report problems to Bruce Ravel (bravel AT bnl DOT gov)

Patches are welcome.

# AUTHOR

Bruce Ravel (bravel AT bnl DOT gov)

[https://github.com/bruceravel](https://github.com/bruceravel)

# LICENSE AND COPYRIGHT

To the extent possible, the authors have waived all rights granted by
copyright law and related laws for the code and documentation that
make up the Perl Interface to the feffpath library.  While information
about Authorship may be retained in some files for historical reasons,
this work is hereby placed in the Public Domain.  This work is
published from: United States.

Note that the feffpath and onepath libraries themselves are NOT public
domain, nor is the Fortran source code for Feff that it relies upon.

Author: Bruce Ravel (bravel AT bnl DOT gov).
Last update: 13 February, 2015
