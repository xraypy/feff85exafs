package Xray::FeffPath;

use Moose;
use Xray::FeffPathWrapper;
use List::MoreUtils qw(any);

has 'wrapper' => (is => 'ro', isa => 'Xray::FeffPathWrapper', default => sub{ Xray::FeffPathWrapper::FEFFPATH->new() });

## simple scalars
has 'Index'   => (is => 'rw', isa => 'Int',  default => 0, trigger => sub{pushback(@_, 'Index'  )},);
has 'nleg'    => (is => 'rw', isa => 'Int',  default => 0, trigger => sub{pushback(@_, 'nleg'   )},);
has 'deg'     => (is => 'rw', isa => 'Num',  default => 0, trigger => sub{pushback(@_, 'deg'    )},);
has 'iorder'  => (is => 'rw', isa => 'Int',  default => 0, trigger => sub{pushback(@_, 'iorder' )},);
has 'nnnn'    => (is => 'rw', isa => 'Bool', default => 0, trigger => sub{pushback(@_, 'nnnn'   )},);
has 'json'    => (is => 'rw', isa => 'Bool', default => 0, trigger => sub{pushback(@_, 'json'   )},);
has 'verbose' => (is => 'rw', isa => 'Bool', default => 0, trigger => sub{pushback(@_, 'verbose')},);
has 'ipol'    => (is => 'rw', isa => 'Bool', default => 0, trigger => sub{pushback(@_, 'ipol'   )},);
has 'elpty'   => (is => 'rw', isa => 'Num',  default => 0, trigger => sub{pushback(@_, 'elpty'  )},);
has 'ne'      => (is => 'rw', isa => 'Int',  default => 0, trigger => sub{pushback(@_, 'ne'     )},);

## arrays
has 'evec'    => (traits  => ['Array'],
		  is      => 'rw',
		  isa     => 'ArrayRef[Str]',
		  default => sub { [0,0,0] },
		  trigger => \&evec_set,
		 );
has 'xivec'   => (traits  => ['Array'],
		  is      => 'rw',
		  isa     => 'ArrayRef[Str]',
		  default => sub { [0,0,0] },
		  trigger => \&xivec_set);

has 'ri'       => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_ri => 'clear', push_ri  => 'push', });
has 'beta'     => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_beta => 'clear', push_beta  => 'push', });
has 'eta'      => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_eta => 'clear', push_eta  => 'push', });

has 'k'        => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_k => 'clear', push_k  => 'push', });
has 'real_phc' => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_real_phc => 'clear', push_real_phc  => 'push', });
has 'mag_feff' => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_mag_feff => 'clear', push_mag_feff  => 'push', });
has 'pha_feff' => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_pha_feff => 'clear', push_pha_feff  => 'push', });
has 'red_fact' => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_red_fact => 'clear', push_red_fact  => 'push', });
has 'lam'      => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_lam => 'clear', push_lam  => 'push', });
has 'rep'      => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
		   handles => {clear_rep => 'clear', push_rep  => 'push', });

sub BUILD {
  my ($self) = @_;
  $self->create_path;
  return $self;
};

## this is the trigger for all the scalar-valued attributes.  it pushes the value back to the wrapper object.
sub pushback {
  my ($self, $new, $old, $which) = @_;
  return if (any {$_ eq $which} qw(nleg ne));
  my $method = join("_", "swig", lc($which), "set");
  $self->wrapper->$method($new);
  $self->ipol(1) if (($which eq 'elpty') and $self->wrapper->swig_elpty_get);
};

## these set the 3-element vectors for polarization and ellipticity
sub evec_set {
  my ($self, $new, $old) = @_;
  $self->wrapper->set_evec(0, $new->[0]);
  $self->wrapper->set_evec(1, $new->[1]);
  $self->wrapper->set_evec(2, $new->[2]);
  $self->ipol($new->[0] or $new->[1] or $new->[2]);
  return $self;
};
sub xivec_set {
  my ($self, $new, $old) = @_;
  $self->wrapper->set_xivec(0, $new->[0]);
  $self->wrapper->set_xivec(1, $new->[1]);
  $self->wrapper->set_xivec(2, $new->[2]);
  $self->ipol($new->[0] or $new->[1] or $new->[2]);
  return $self;
};


sub create_path {
  my ($self) = @_;
  $self->wrapper->create_path;
  foreach my $att (qw(Index nleg deg iorder nnnn json verbose ipol elpty)) {
    my $method = join("_", "swig", lc($att), "get");
    $self->$att($self->wrapper->$method||0);
  };
  return $self;
};

sub atom {
  my ($self, $x, $y, $z, $ip) = @_;
  my $nleg = $self->wrapper->add_scatterer($x, $y, $z, $ip);
  $self->nleg($nleg);
  return $self;
};

sub path {
  my ($self) = @_;
  $self->wrapper->make_path;
  ## clear out all the arrays
  $self->ne($self->wrapper->swig_ne_get);
  $self->clear_ri;
  $self->clear_beta;
  $self->clear_eta;
  $self->clear_k;
  $self->clear_real_phc;
  $self->clear_mag_feff;
  $self->clear_pha_feff;
  $self->clear_red_fact;
  $self->clear_rep;
  $self->clear_lam;
  ## fill then with the results of the calculation
  for my $i (0 .. $self->nleg-1) {
    $self->push_ri($self->wrapper->get_ri($i));
    $self->push_beta($self->wrapper->get_beta($i));
    $self->push_eta($self->wrapper->get_eta($i));
  };
  for my $i (0 .. $self->ne-1) {
    $self->push_k($self->wrapper->get_k($i));
    $self->push_real_phc($self->wrapper->get_real_phc($i));
    $self->push_mag_feff($self->wrapper->get_mag_feff($i));
    $self->push_pha_feff($self->wrapper->get_pha_feff($i));
    $self->push_red_fact($self->wrapper->get_red_fact($i));
    $self->push_lam($self->wrapper->get_lam($i));
    $self->push_rep($self->wrapper->get_rep($i));
  };
  return $self;
};

sub clear {
  my ($self) = @_;
  $self->wrapper->clear_path;
  foreach my $att (qw(Index nleg deg iorder nnnn json verbose ipol elpty)) {
    my $method = join("_", "swig", lc($att), "get");
    $self->$att($self->wrapper->$method||0);
  };
  $self->evec([0,0,0]);
  $self->xivec([0,0,0]);
  $self->clear_ri;
  $self->clear_beta;
  $self->clear_eta;
  $self->clear_k;
  $self->clear_real_phc;
  $self->clear_mag_feff;
  $self->clear_pha_feff;
  $self->clear_red_fact;
  $self->clear_rep;
  $self->clear_lam;
};

# use Term::ANSIColor qw(:constants);
# sub trace {
#   my ($self) = @_;
#   my $max_depth = 30;
#   my $i = 0;
#   my $base = substr($INC{'Demeter.pm'}, 0, -10);
#   my ($green, $red, $yellow, $end) = (BOLD.GREEN, BOLD.RED, BOLD.YELLOW, RESET);
#   local $|=1;
#   print($/.BOLD."--- Begin stack trace ---$end\n");
#   while ( (my @call_details = (caller($i++))) && ($i<$max_depth) ) {
#     (my $from = $call_details[1]) =~ s{$base}{};
#     my $line  = $call_details[2];
#     my $color = RESET.YELLOW;
#     (my $func = $call_details[3]) =~ s{(?<=::)(\w+)\z}{$color$1};
#     print("$green$from$end line $red$line$end in function $yellow$func$end\n");
#   }
#   print(BOLD."--- End stack trace ---$end\n");
#   return $self;
# };



no Moose;
__PACKAGE__->meta->make_immutable;


=head1 NAME

Xray::FeffPath - SWIG-based interface to libfeffpath

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
  use Xray::FeffPath;

  my $path = Xray::FeffPath->new();
  $path->create_path;
  $path->atom(0, 0, -3.61, 1);
  $path->atom(-1.805, 0, -1.805, 1);
  $path->path;
  my @kgrid  = $path->k;
  my @caps   = $path->real_phc;
  my @amff   = $path->mag_feff;
  my @phff   = $path->pha_feff;
  my @redfac = $path->red_fact;
  my @lambda = $path->lam;
  my @rep    = $path->rep;
  undef $path

=head1 METHODS

=over 4

=item C<new>

Create the FeffPath object

   my $path = Xray::FeffPath->new();

=item C<create_path>

Initialize the underlying C struct.  This is the necessary first step
to a feffpath calculation.

   $path->create_path;

=item C<atom>

Add an atom to the scattering path, supplying the Cartesian
coordiantes of the atom and its unique potential index.

   $nleg = $path->atom($x, $y, $z, $ipot);

This returns the current number of legs in the path.  To make a SS
paths do:

   $path->atom($x,  $y,  $z, $ipot1);
   print $path->nleg, $/;
   #   ==prints===> 2

To make a MS paths do:

   $path->atom($x1, $y1, $z1, $ipot1);
   $path->atom($x2, $y2, $z2, $ipot2);
   $path->atom($x3, $y3, $z3, $ipot3);
   print $path->nleg, $/;
   #   ==prints===> 4

Note that there is no setter method for the nleg attribute.  It is set
by the call to the C<atom> method.

=item C<path>

Compute the scattering path

   $path->path;

=item C<clear>

Reinitialize the scattering path

   $path->clear;

=back

=head1 ATTRIBUTES

Each of the members of the FEFFPATH struct is wrapped in a getter
method of the same name.  Each of the following exists:

=head2 Scalar valued attributes

=over 4

=item C<Index> (integer, default = 9999)

The path index, used for form the name of the output F<feffNNNN.dat>
file.

Accessing this and other scalar-values attributes works in the typical
Moose-y way:

   $path->Index(4);
   print $path->Index, $/;
   ## ==prints==> 4

=item C<deg> (float, default = 0)

The path degeneracy.  This is required input for a calculation.

=item C<nleg> (integer, default = 0)

The number of legs in the path.  This should never be set by hand.  It
gets set correctly whenever the C<atom> method is called.

=item C<iorder> (integer, default = 2)

Order of approximation for multiple scattering paths.

=item C<nnnn> (boolean, default = 0)

Flag for writing F<feffNNNN.dat> file.

=item C<json> (boolean, default = 0)

Flag for writing F<feffNNNN.json> file.

=item C<verbose> (integer, default = 0)

Flag for writing screen messages as files are written.

=item C<ipol> (boolean, default = 0)

Perform polarization calculation.

=item C<elpty> (integer, default = 0)

Ellipticity for use in polarization calculation.

=item C<reff> (float, default = 0)

Half path length.  This should never be set by hand.  It
gets set correctly whenever the C<path> method is called.

=item C<ne> (integer, default = 0)

Number of outpit k grid point.  This should never be set by hand.  It
gets set correctly whenever the C<path> method is called.

=back

=head2 Array reference valued attributes

=over 4

=item C<evec>

The electric vector of the incident beam for the polarization
calculation.  Note that the setter's argument is an array reference as
is the getter's return value.

    $path->evec([0,0,1]);
    $evec_ref = $path->evec;
    print join(", ", @$evec_ref);
    ## ==prints==> 0, 0, 1

=item C<xivec>

The pointing vector of the incident beam for the ellipticity
calculation.  Note that the setter's argument is an array reference as
is the getter's return value.

    $path->xivec([1,1,0]);
    $xivec_ref = $path->xivec;
    print join(", ", @$xivec_ref);
    ## ==prints==> 1, 1, 0

=item C<ri>

The leg lengths (in Angstroms) of the specified path.  This gets set
when C<path> is called.  It returns an reference to an array with
C<nleg> elements.

=item C<beta>

The Eulerian beta angles (in degrees) of the specified path.  This
gets set when C<path> is called.  It returns an reference to an array
with C<nleg> elements.

=item C<eta>

The Eulerian eta angles (in degrees) of the specified path.  This gets
set when C<path> is called.  It returns an reference to an array with
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

A reference to an array containing lambda, themean free path.  This is
column 6 from F<feffNNNN.dat>.

=item C<rep>

A reference to an array containing the real part of the comple
momentum.  This is column 7 from F<feffNNNN.dat>.


=back

=head1 DEPENDENCIES

=over 4

=item *

L<Moose>

=back

=head1 BUGS AND LIMITATIONS

=over 4

=item *

Getter method(s) for rat and ipot

=item *

Error handling is almost nonexistent

=back

Please report problems to Bruce Ravel (bravel AT bnl DOT gov)

Patches are welcome.

=head1 AUTHOR

Bruce Ravel (bravel AT bnl DOT gov)

L<http://cars9.uchicago.edu/~ravel/software/>

=head1 LICENSE AND COPYRIGHT

To the extent possible, the authors have waived all rights granted by
copyright law and related laws for the code and documentation that
make up the Perl Interface to the feffpath library.  While information
about Authorship may be retained in some files for historical reasons,
this work is hereby placed in the Public Domain.  This work is
published from: United States.

Note that the feffpath library itself is NOT public domain, nor is the
Fortran source code it relies upon.

Author: Bruce Ravel (bravel AT bnl DOT gov).
Last update: 2 October, 2014

=cut
