package Xray::Feff::Phases;

use Moose;
use MooseX::NonMoose;
use MooseX::Aliases;
extends 'Xray::Feff::PhasesWrapper';

use List::MoreUtils qw(any);

our $VERSION = '1.00'; # Inline::MakeMake uses /^\d.\d\d$/ as the pattern for the version number -- note the two digits to the right of the dot

has 'wrapper' => (
		  is        => 'ro',
		  #traits => [qw(NoClone)],
		  isa       => 'Xray::Feff::PhasesWrapper',
		  init_arg  => undef,
		  default   => sub{ Xray::Feff::PhasesWrapper->new() },
		  #lazy      => 1,
		  #builder   => '_build_object',
		 );

has 'errorcode'    => (is => 'rw', isa => 'Int',  default => 0,);
has 'errormessage' => (is => 'rw', isa => 'Str',  default => q{},);

## constants from feffpath.h
has 'natx'     => (is => 'ro', isa => 'Int', default => sub{ Xray::Feff::PhasesWrapper->_natx   });
has 'nphx'     => (is => 'ro', isa => 'Int', default => sub{ Xray::Feff::PhasesWrapper->_nphx   });
has 'nheadx'   => (is => 'ro', isa => 'Int', default => sub{ Xray::Feff::PhasesWrapper->_nheadx });
has 'bohr'     => (is => 'ro', isa => 'Num', default => sub{ Xray::Feff::PhasesWrapper->_bohr   });
has 'ryd'      => (is => 'ro', isa => 'Num', default => sub{ Xray::Feff::PhasesWrapper->_ryd    });
has 'hart'     => (is => 'ro', isa => 'Num', default => sub{ Xray::Feff::PhasesWrapper->_hart   });


## scalars
has 'jsonfile' => (is => 'rw', isa => 'Str',  default => 'libpotph.json', trigger => sub{pushback(@_, 'jsonfile' )});
has 'ntitle'   => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'ntitle'   )});
has 'nat'      => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'nat'      )});
has 'nph'      => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'nph'      )});
has 'ihole'    => (is => 'rw', isa => 'Int',  default => 1,               trigger => sub{pushback(@_, 'ihole'    )});
has 'rscf'     => (is => 'rw', isa => 'Num',  default => 4.0,             trigger => sub{pushback(@_, 'rscf'     )});
has 'lscf'     => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'lscf'     )});
has 'nscmt'    => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'nscmt'    )});
has 'ca'       => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'ca'       )});
has 'nmix'     => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'nmix'     )});
has 'ecv'      => (is => 'rw', isa => 'Num',  default => -40.0,           trigger => sub{pushback(@_, 'ecv'      )});
has 'icoul'    => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'icoul'    )});
has 'ipol'     => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'ipol'     )});
has 'elpty'    => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'elpty'    )});
has 'ispin'    => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'ispin'    )});
has 'angks'    => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'angks'    )});
has 'gamach'   => (is => 'rw', isa => 'Num',  default => 0.0);
has 'ixc'      => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'ixc'      )});
has 'vr0'      => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'vr0'      )});
has 'vi0'      => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'vi0'      )});
has 'ixc0'     => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'ixc0'     )});
has 'iafolp'   => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'iafolp'   )});
has 'rgrd'     => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'rgrd'     )});
has 'iunf'     => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'iunf'     )});
has 'inters'   => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'inters'   )});
has 'totvol'   => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'totvol'   )});
has 'jumprm'   => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'jumprm'   )});
has 'nohole'   => (is => 'rw', isa => 'Int',  default => 0,               trigger => sub{pushback(@_, 'nohole'   )});


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
has 'spvec'    => (traits  => ['Array'],
		   is      => 'rw',
		   isa     => 'ArrayRef[Num]',
		   default => sub { [0,0,0] },
		   trigger => \&spvec_set, );


has 'iz'       => (is => 'rw', isa => 'ArrayRef', default => sub{[]}); # trigger => \&iz_set
has 'lmaxsc'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'lmaxph'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'xnatph'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'spinph'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'folp'     => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'xion'     => (is => 'rw', isa => 'ArrayRef', default => sub{[]});


has 'iphat'    => (is => 'rw', isa => 'ArrayRef', default => sub{[]});


# has 'rat'      => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
# has 'potlbl'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
# has 'titles'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]});



sub BUILD {
  my ($self) = @_;
  $self->wrapper->_create_phases;
  return $self;
};

sub DEMOLISH {
  my ($self) = @_;
  #print "in my DEMOLISH\n";
  #print $self->wrapper, $/;
  $self->wrapper->_cleanup;
  return $self;
};


sub _json {
  my ($self) = @_;
  my $ret = $self->wrapper->_read_libpotph_json;
  $self->_fetch;
  return $ret;
};

sub phases {
  my ($self) = @_;
  my $ret = $self->wrapper->_make_phases;
  return $ret;
};

sub _fetch {
  my ($self) = @_;
  foreach my $att (qw(ntitle nat nph ihole rscf lscf nscmt ca nmix ecv icoul elpty ispin angks gamach
		      ixc vr0 vi0 ixc0 iafolp rgrd iunf inters totvol jumprm nohole)) {
    my $method = '_'.$att;
    $self->$att($self->wrapper->$method);
  };
  print $self->wrapper->_iz_array, $/;
  return $self;
};


## ------------------------------------------------------------
## trigger methods

sub pushback {
  my ($self, $new, $old, $which) = @_;
  return if (any {$_ eq $which} qw(gamach));
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
  $self->ipol(1);
  return $self;
};
sub xivec_set {
  my ($self, $vec) = @_;
  $self->wrapper->_set_xivec(1.0*$vec->[0], 1.0*$vec->[1], 1.0*$vec->[2]);
  $self->ipol(1);
  return $self;
};
sub spvec_set {
  my ($self, $vec) = @_;
  $self->wrapper->_set_spvec(1.0*$vec->[0], 1.0*$vec->[1], 1.0*$vec->[2]);
  $self->ispin(1);
  return $self;
};







no Moose;
__PACKAGE__->meta->make_immutable;
1;



=head1 NAME

Xray::FeffPath - Moose wrapper around an Inline::C interface to libfeffpath

=head1 VERSION

1.00

=head1 SYNOPSIS

This provides a simple object oriented interface to the C<feffphases>
wrapper around Feff's C<libpotph> subroutine.  C<feffphases> combines
the potentials and phase shifts calculations, writing a F<phase.pad>
file to disk.

The following computes the phase shifts for copper metal:

  #!/usr/bin/perl
  use strict;
  use warnings;
  use Xray::Feff::Phases;

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

Create the FeffPath object

   my $phases = Xray::Feff::Phases->new();

=item C<clear>

Reinitialize the scattering path

   $path->clear;

=back

To destroy a FeffPath object, simply

   undef $path;


=head1 ATTRIBUTES

Each of the members of the FEFFPATH struct is wrapped in a getter
method of the same name.  Each of the following exists:

=head2 Scalar valued attributes

=over 4

=item C<jsonfile> (character, default = libpotph.json)

The path to the F<libpth.json> file.




=item C<errorcode> (integer)

An integer error code from either C<phases>.  Check for non-zero
value.  See C<errormessage> for an explanation of the problem.

=item C<errormessage> (string)

A explanation of words of the problem found during C<phases>.

   $path->phases;
   if ($path->errorcode) {
      print $path->errormessage;
   };
   ## ==prints==>
    Blah blah blah

=back


=head2 Constants from feffpath.h

The following constants are captured from the feffpath header file:

=over 4

=item C<natx> (1000)

The maximum number of atoms in the cluster.

=item C<nphx> (11)

The maximum number of unique potentials.

=item C<nheadx> (30)

The maximum number of title l;ines

=item C<bohr> (0.529177249)

A unit of length, in Angstrom.

=item C<ryd> (13.605698)

A Rydberg, a unit of energy, in eV.

=item C<hart> (27.211396)

A Hartree, a unit of energy, in eV.

=back

=head2 Attributes with values of array-reference

=over 4

=item C<evec>

The electric vector of the incident beam for the polarization
calculation.  Note that the setter's argument is an array reference as
is the getter's return value.

    $phases->evec([0,0,1]);
    $evec_ref = $phases->evec;
    print join(", ", @$evec_ref);
    ## ==prints==> 0, 0, 1

Setting C<evec> will also set C<ipol> to true.

=item C<xivec>

The Poynting vector of the incident beam for the ellipticity
calculation.  Note that the setter's argument is an array reference as
is the getter's return value.

    $phases->xivec([1,1,0]);
    $xivec_ref = $phases->xivec;
    print join(", ", @$xivec_ref);
    ## ==prints==> 1, 1, 0

Setting C<xivec> will also set C<ipol> to true.

=item C<spvec>

The Poynting vector of the incident beam for the ellipticity
calculation.  Note that the setter's argument is an array reference as
is the getter's return value.

    $phases->xivec([1,1,0]);
    $xivec_ref = $phases->xivec;
    print join(", ", @$xivec_ref);
    ## ==prints==> 1, 1, 0

Setting C<xivec> will also set C<ispin> to true.


=item C<iz>

A reference to an array containing the Z numbers of each unique
potential.

=back

=head1 Error codes

The error codes returned by the C<phases> methods are the sums of the
codes actually found by the method call, thus the returned error codes
are meant to be interpreted bitwise.

=over 4

=item I<1>

=item I<2>

=item I<4>

=item I<8>

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

Blah blah

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

Note that the feffphases and libpotph libraries themselves are NOT public
domain, nor is the Fortran source code for Feff that it relies upon.

Author: Bruce Ravel (bravel AT bnl DOT gov).
Created: 5 June, 2015

=cut
