package Xray::Feff::Path;

use Moose;
use MooseX::NonMoose;
use MooseX::Aliases;
extends 'Xray::Feff::PathWrapper';

use List::MoreUtils qw(any);

our $VERSION = '1.00'; # Inline::MakeMake uses /^\d.\d\d$/ as the
                       # pattern for the version number -- note the
                       # two digits to the right of the dot

has 'wrapper' => (
		  is        => 'ro',
		  #traits => [qw(NoClone)],
		  isa       => 'Xray::Feff::PathWrapper',
		  init_arg  => undef,
		  default   => sub{ Xray::Feff::PathWrapper->new() },
		  #lazy      => 1,
		  #builder   => '_build_object',
		 );

has 'errorcode'    => (is => 'rw', isa => 'Int',  default => 0,);
has 'errormessage' => (is => 'rw', isa => 'Str',  default => q{},);

## constants from feffpath.h
has 'nex'      => (is => 'ro', isa => 'Int', default => sub{ Xray::Feff::PathWrapper->_nex    });
has 'nphx'     => (is => 'ro', isa => 'Int', default => sub{ Xray::Feff::PathWrapper->_nphx   });
has 'npatx'    => (is => 'ro', isa => 'Int', default => sub{ Xray::Feff::PathWrapper->_npatx  });
has 'legtot'   => (is => 'ro', isa => 'Int', default => sub{ Xray::Feff::PathWrapper->_legtot });
has 'bohr'     => (is => 'ro', isa => 'Num', default => sub{ Xray::Feff::PathWrapper->_bohr   });
has 'ryd'      => (is => 'ro', isa => 'Num', default => sub{ Xray::Feff::PathWrapper->_ryd    });
has 'hart'     => (is => 'ro', isa => 'Num', default => sub{ Xray::Feff::PathWrapper->_hart   });

## location of phase.pad file
has 'phpad'    => (is => 'rw', isa => 'Str',  default => 'phase.pad', trigger => sub{pushback(@_, 'phpad'  )},
		   alias => 'phbin');

## basic path parameters
has 'Index'    => (is => 'rw', isa => 'Int',  default => 9999, trigger => sub{pushback(@_, 'Index'  )},);
has 'nleg'     => (is => 'rw', isa => 'Int',  default => 0,    trigger => sub{pushback(@_, 'nleg'   )},);
has 'degen'    => (is => 'rw', isa => 'Num',  default => 1.0,  trigger => sub{pushback(@_, 'degen'  )},
		   alias => 'degeneracy');
has 'iorder'   => (is => 'rw', isa => 'Int',  default => 2,    trigger => sub{pushback(@_, 'iorder' )},);

## output and verbosity control
has 'nnnn'     => (is => 'rw', isa => 'Bool', default => 0,    trigger => sub{pushback(@_, 'nnnn'   )},);
has 'xdi'      => (is => 'rw', isa => 'Bool', default => 0,    trigger => sub{pushback(@_, 'xdi'    )},);
has 'verbose'  => (is => 'rw', isa => 'Bool', default => 0,    trigger => sub{pushback(@_, 'verbose')},);

## header information from feffNNNN.dat
has 'edge'     => (is => 'rw', isa => 'Num',  default => 0.0);
has 'gam_ch'   => (is => 'rw', isa => 'Num',  default => 0.0);
has 'kf'       => (is => 'rw', isa => 'Num',  default => 0.0);
has 'mu'       => (is => 'rw', isa => 'Num',  default => 0.0);
has 'rnorman'  => (is => 'rw', isa => 'Num',  default => 0.0);
has 'rs_int'   => (is => 'rw', isa => 'Num',  default => 0.0);
has 'vint'     => (is => 'rw', isa => 'Num',  default => 0.0);
has 'exch'     => (is => 'rw', isa => 'Str',  default => '');
has 'version'  => (is => 'rw', isa => 'Str',  default => '');

## polarization
has 'ipol'     => (is => 'rw', isa => 'Int',  default => 0,    trigger => sub{pushback(@_, 'ipol'   )},);
has 'elpty'    => (is => 'rw', isa => 'Num',  default => 0.0,  trigger => sub{pushback(@_, 'elpty'  )},
		   alias => 'ellipticity');
has 'evec'     => (traits  => ['Array'],
		   is      => 'rw',
		   isa     => 'ArrayRef[Num]',
		   default => sub { [0,0,0] },
		   trigger => \&evec_set, );
has 'xivec'    => (traits  => ['Array'],
		   is      => 'rw',
		   isa     => 'ArrayRef[Num]',
		   default => sub { [0,0,0] },
		   trigger => \&xivec_set, );

## geometry table
has 'absorber' => (is => 'rw', isa => 'ArrayRef', default => sub{[0,0,0]});
has 'ipot'     => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Int]', default => sub { [] },
                   handles => {clear_ipot => 'clear', push_ipot  => 'push', });
has 'rat'      => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[ArrayRef[Num]]', default => sub { [] },
                   handles => {clear_rat => 'clear', push_rat  => 'push', });
has 'iz'       => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Int]', default => sub { [] },
                   handles => {clear_iz => 'clear', push_iz  => 'push', });
has 'reff'     => (is => 'rw', isa => 'Num',  default => 0.0);
has 'ri'       => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
                   handles => {clear_ri => 'clear', push_ri => 'push', });
has 'beta'     => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
                   handles => {clear_beta => 'clear', push_beta => 'push', });
has 'eta'      => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
                   handles => {clear_eta => 'clear', push_eta  => 'push', });

## data table
has 'ne'       => (is => 'rw', isa => 'Int',      default => 0);
has 'k'        => (is => 'rw', isa => 'ArrayRef', default => sub{[]}, alias => 'kgrid');
has 'real_phc' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'mag_feff' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'pha_feff' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'red_fact' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'lam'      => (is => 'rw', isa => 'ArrayRef', default => sub{[]}, alias => 'lambda');
has 'rep'      => (is => 'rw', isa => 'ArrayRef', default => sub{[]}, alias => 'realp');


sub BUILD {
  my ($self) = @_;
  $self->wrapper->_create_path;
  return $self;
};

sub DEMOLISH {
  my ($self) = @_;
  #print "in my DEMOLISH\n";
  #print $self->wrapper, $/;
  $self->wrapper->_cleanup;
  return $self;
};


## ------------------------------------------------------------
## trigger methods

sub pushback {
  my ($self, $new, $old, $which) = @_;
  return if (any {$_ eq $which} qw(nleg ne edge gam_ch kf mu rnorman version exch rs_int vint));
  my $method = '_set_' . lc($which);
  if ($self->meta->get_attribute($which)->type_constraint eq 'Num') {
    $self->wrapper->$method(1.0*$new);
  } elsif ($self->meta->get_attribute($which)->type_constraint eq 'Bool') {
    my $val = ($new) ? 1 : 0;
    $self->wrapper->$method($val);
  } else {
    $self->wrapper->$method($new);
  };
  $self->ipol(1) if (($which eq 'elpty') and $self->wrapper->_elpty);
};

sub evec_set {
  my ($self, $vec) = @_;
  $self->wrapper->_set_evec(1.0*$vec->[0], 1.0*$vec->[1], 1.0*$vec->[2]);
  return $self;
};
sub xivec_set {
  my ($self, $vec) = @_;
  $self->wrapper->_set_xivec(1.0*$vec->[0], 1.0*$vec->[1], 1.0*$vec->[2]);
  return $self;
};


## ------------------------------------------------------------
## main methods

sub clear {
  my ($self) = @_;
  $self->wrapper->_clear_path;
  foreach my $att (map {$_->name} $self->meta->get_all_attributes) {              # iterate over all attribute names
    next if (not $self->meta->get_attribute($att)->get_write_method);             # skip ro attributes, see Class/MOP/Attribute.pm#Informational
    next if ($self->meta->get_attribute($att)->type_constraint =~ m{Feff::Path});   # skip wrapper
    my $method = '_'.lc($att);

    if ($att =~ m{(?:(?:e|xi)vec|absorber)}) {					  # 3-vec attributes
      $self->$att([0,0,0])

    } elsif ($self->meta->get_attribute($att)->type_constraint =~ m{ArrayRef}) {  # array valued attributes
      $self->$att([])

    } elsif ($self->meta->get_attribute($att)->type_constraint =~ m{Bool}) {	  # boolean valued attributes
      my $val = $self->wrapper->$method || 0;
      $self->$att($val);

    } else {									  # all the rest
      my $val = $self->wrapper->$method || 0;
      $self->$att($val);
    };
  };
  return $self;
};

sub atom {
  my ($self, $x, $y, $z, $ip) = @_;
  my $message;
  my $xx  = $x - $self->absorber->[0];
  my $yy  = $y - $self->absorber->[1];
  my $zz  = $z - $self->absorber->[2];
  my $err = $self->wrapper->_add_scatterer($xx, $yy, $zz, $ip);
  $self->errorcode($err);
  my $em  = $self->wrapper->_errormessage;
  $em =~ s{add_scatterer}{atom method};
  $self->errormessage($em);
  $self->nleg($self->wrapper->_nleg);
  $self->push_ipot($ip);
  $self->push_rat([$x, $y, $z]);
  return $err;
};

sub path {
  my ($self) = @_;
  $self->push_ipot(0);
  $self->push_rat($self->absorber);
  my $err = $self->wrapper->_make_path;
  $self->errorcode($err);
  my $em = $self->wrapper->_errormessage;
  $em =~ s{make_path}{path method};
  $self->errormessage($em);
  if (not $err) {

    ## grab scalar output
    foreach my $a (qw(reff ne edge gam_ch kf mu rnorman rs_int vint exch version)) {
      my $method = '_'.$a;
      $self->$a($self->wrapper->$method);
    };

    ## grab data table
    foreach my $a (qw(k real_phc mag_feff pha_feff red_fact lam rep)) {
      my $method = '_' . $a . '_array';
      my @array = $self->wrapper->$method;
      $self->$a(\@array);
    };

    ## grab geometry table
    foreach my $a (qw(ri beta eta iz)) {
      my $method = '_' . $a . '_array';
      my @array = $self->wrapper->$method;
      $self->$a(\@array);
    };
  };
  return $err;
};


no Moose;
__PACKAGE__->meta->make_immutable;
1;



=head1 NAME

Xray::Feff::Path - Moose wrapper around an Inline::C interface to libfeffpath

=head1 VERSION

1.00

=head1 SYNOPSIS

This provides a simple object oriented interface to the C<feffpath>
wrapper around Feff's C<onepath> subroutine.  C<onepath> combines the
F-matrix calculation in F<GENFMT/genfmt.f> with the construction of the
common presentation of F-effective from F<FF2X/feffdt.f>.

Given the geometry of a scattering path in the form of Cartesian
coordinates of atoms making up the path, along with some additional
information such as the degeneracy, compute F-effective for that
scattering path.

This replaces the need to read a F<feffNNNN.dat> file.  Given a
scattering geometry, the equivalent data can be computed directly.

The following computes the columns of the F<feff0004.dat> file for
copper metal:

  #!/usr/bin/perl
  use strict;
  use warnings;
  use Xray::Feff::Path;

  my $path = Xray::Feff::Path->new();
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

=head1 INSTALLATION

After you have built and installed I<feff85exafs>, do the following:

  perl Makefile.PL
  make
  make test
  sudo make install

That's it!  Note, though, that building this wrapper B<requires> that
the fortran and C code I<feff85exafs> be completely compiled and that
the resulting libraries (and other files) be successfully installed.

=head1 METHODS

=over 4

=item C<new>

Create the Feff::Path object

   my $path = Xray::Feff::Path->new();

=item C<atom>

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
C<nleg> is incremented by each call to the C<atom> method.

=item C<path>

Compute the scattering path

   $path->path;

=item C<clear>

Reinitialize the scattering path

   $path->clear;

=back

To destroy a Feff::Path object, simply

   undef $path;

=head1 ATTRIBUTES

Each of the members of the FEFFPATH struct is wrapped in a getter
method of the same name.  Each of the following exists:

=head2 Scalar valued attributes

=over 4

=item C<phpad> (character, default = phase.pad)

The path to the F<phase.pad> file.  This gets right-padded with spaces
to 256 characters, which is the length of the parameter in the
underlying Fortran library.  Therefore, the value returned by the
getter will be right-padded.

=item C<Index> (integer, default = 9999)

The path index, used for form the name of the output F<feffNNNN.dat>
file.

Accessing this and other scalar-values attributes works in the typical
Moose-y way:

   $path->Index(4);
   print $path->Index, $/;
   ## ==prints==> 4

Note that this attribute (but only this one) is capitalized to avoid
confusion with the built-in C<index> function.

=item C<degen> (float, default = 0)

The path degeneracy.  This is required input for a calculation.

=item C<nleg> (integer, default = 0)

The number of legs in the path.  This should never be set by hand.  It
gets incremented by calls to the C<atom> method.

=item C<iorder> (integer, default = 2)

Order of approximation for multiple scattering paths.

=item C<nnnn> (boolean, default = 0)

Flag for writing F<feffNNNN.dat> file.

=item C<xdi> (boolean, default = 0)

Flag for writing F<feffNNNN.xdi> file.

=item C<verbose> (integer, default = 0)

Flag for writing screen messages as files are written.

=item C<ipol> (boolean, default = 0)

Perform polarization calculation.  This gets toggled to true whenever
any of C<elpty>, C<evec>, or C<xivec> are set.

=item C<elpty> (integer, default = 0)

Ellipticity for use in polarization calculation.

=item C<reff> (float, default = 0)

Half path length.  This should never be set by hand.  It
gets set correctly when the C<path> method is called.

=item C<ne> (integer, default = 0)

Number of outpit k grid point.  This should never be set by hand.  It
gets set correctly when the C<path> method is called.

=back

=head2 Error reporting

=over 4

=item C<errorcode> (integer)

An integer error code from either C<atom> or C<path>.  Check for
non-zero value.  For a non-zero return from C<atom>, calling C<path>
is likely to cause a problem in the genfmt calculation.  For a
non-zero return from C<path>, the genfmt calculation was skipped.  See
C<errormessage> for an explanation of the problem.

=item C<errormessage> (string)

A explanation of words of the problem found during C<atom> or C<path>.

   $path->degen(-48);
   $path->atom(0, 0, 3.61, 1);
   $path->path;
   if ($path->errorcode) {
      print $path->errormessage;
   };
   ## ==prints==>
    Error in make_path
        path degeneracy (-48.00) is negative

=back

=head2 Attributes related to the potential model

Several parameters are captured by the C struct (and therefore
available via the perl wrapper) which are related to how Feff's
potential model was calculated.  This information is written to the
F<feffNNNN.dat> header.

=over 4

=item  C<version>

A string identifying the version number and the `feffpath` wrapper
revision.

=item  C<exch>

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
L<http://feffproject.org/feff/Docs/feff9/feff90/feff90_users_guide.pdf>
at page 92.

=item  C<edge>

A poor estimate of the energy threshold relative to the atomic value, in eV.

=item  C<gam_ch>

The core-hole lifetime, in eV.

=item  C<kf>

The k value at the Fermi energy, in inverse Angstroms.

=item  C<mu>

The Fermi energy, in eV.

=item  C<rnorman>

The average Norman radius.

=item  C<rs_int>

An estimate of the interstitial radius.

=item  C<vint>

The interstitial potential, in eV.

=back

=head2 Constants from feffpath.h

The following constants are captured from the feffpath header file:

=over 4

=item C<nex> (150)

The maximum number of points in the energy grid.

=item C<nphx> (11)

The maximum number of unique potentials.

=item C<npatx> (8)

The maximum number of path atoms

=item C<legtot> (9)

The maximum number of legs in a path (npatx + 1)

=item C<bohr> (0.529177249)

A unit of length, in Angstrom.

=item C<ryd> (13.605698)

A Rydberg, a unit of energy, in eV.

=item C<hart> (27.211396)

A Hartree, a unit of energy, in eV.

=back

=head2 Attributes with values of array-reference

=over 4

=item C<absorber>

The Cartesian coordinates of the absober atom.  By default the
absorber is placed at (0,0,0).

=item C<ipot>

A list of the unique potentials of the scatterers in a path.  This is
populated by the calls to the C<atom> method.

=item C<rat>

A list of lists of the Cartesian coordinates of the scatterers in a
path.  This is populated by the calls to the C<atom> method.

=item C<evec>

The electric vector of the incident beam for the polarization
calculation.  Note that the setter's argument is an array reference as
is the getter's return value.

    $path->evec([0,0,1]);
    $evec_ref = $path->evec;
    print join(", ", @$evec_ref);
    ## ==prints==> 0, 0, 1

Setting C<evec> will also set C<ipol> to true.

=item C<xivec>

The Poynting vector of the incident beam for the ellipticity
calculation.  Note that the setter's argument is an array reference as
is the getter's return value.

    $path->xivec([1,1,0]);
    $xivec_ref = $path->xivec;
    print join(", ", @$xivec_ref);
    ## ==prints==> 1, 1, 0

Setting C<xivec> will also set C<ipol> to true.

=item C<ri>

The leg lengths (in Angstroms) of the specified path.  This gets set
when C<path> is called.  It returns an reference to an array with
C<nleg> elements.

=item C<beta>

The Eulerian beta angles (in degrees) of the specified path.  This
gets set when C<path> is called.  It returns a reference to an array
with C<nleg> elements.

=item C<eta>

The Eulerian eta angles (in degrees) of the specified path.  This gets
set when C<path> is called.  It returns a reference to an array with
C<nleg> elements.

=item C<k>

A reference to an array containing the k grid of the calculation.
This is column 1 from F<feffNNNN.dat>.

=item C<real_phc>

A reference to an array containing the central atom phase shifts.
This is column 2 from F<feffNNNN.dat>.

=item C<mag_feff>

A reference to an array containing the magnitude of F-effective.  This
is column 3 from F<feffNNNN.dat>.

=item C<pha_feff>

A reference to an array containing the phase of F-effective.  This is
column 4 from F<feffNNNN.dat>.

=item C<red_fact>

A reference to an array containing the reduction factor.  This is
column 5 from F<feffNNNN.dat>.

=item C<lam>

A reference to an array containing lambda, the mean free path.  This is
column 6 from F<feffNNNN.dat>.

=item C<rep>

A reference to an array containing the real part of the complex
momentum.  This is column 7 from F<feffNNNN.dat>.

=back

=head1 Error codes

The error codes returned by the C<atom> and C<path> methods are the
sums of the codes actually found by the method call, thus the returned
error codes are meant to be interpreted bitwise.

=head2 C<atom> method

If any of these are triggered, it is very likely that a call to
C<path> will return unreliable results or may crash the program.

=over 4

=item I<1>

ipot argument to add_scatterer is less than 0

=item I<2>

ipot argument to add_scatterer is greater than 7

=item I<4>

coordinates are for an atom too close to the previous atom in the path

=item I<8>

nlegs greater than legtot

=back

=head2 C<path> method

If any of these are triggered, genfmt will not be called and all
output arrays will be zero-filled.

=over 4

=item I<1>

the first atom specified is the absorber

=item I<2>

the last atom specified is the absorber

=item I<4>

path degeneracy is negative

=item I<8>

path index not between 0 and 9999

=item I<16>

ellipticity not between 0 and 1

=item I<32>

iorder not between 0 and 10

=item I<64>

F<phase.pad> cannot be found or cannot be read

=back


=head1 EXTERNAL DEPENDENCIES

=over 4

=item *

L<Inline>, and L<Inline::C>

=item *

L<Moose>

=item *

L<MooseX::NonMoose> and L<MooseX::Aliases>

=back

=head1 BUGS AND LIMITATIONS

=over 4

=item *

Setting a boolean unsets phpad at Inline::C level.  Weird!
Work-around is at line 130.

=item *

Polarization test is quite trivial.  Should test against an actual
calculation.

=item *

C<ipol> and C<rat> not accessed via wrapper.  That said,
Xray::Feff::Path keeps those attributes current.

=back

Please report problems to Bruce Ravel (bravel AT bnl DOT gov)

Patches are welcome.

=head1 AUTHOR

Bruce Ravel (bravel AT bnl DOT gov)

L<https://github.com/bruceravel>

=head1 LICENSE AND COPYRIGHT

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

=cut
