package Xray::FeffPath;

use Moose;
use MooseX::NonMoose;
use MooseX::Aliases;
extends 'Xray::FeffPathWrapper';

use List::MoreUtils qw(any);

our $VERSION = '1.00'; # Inline::MakeMake uses /^\d.\d\d$/ as the pattern for the version number -- note the two digits to the right of the dot

has 'wrapper' => (
		  is        => 'ro',
		  #traits => [qw(NoClone)],
		  isa       => 'Xray::FeffPathWrapper',
		  init_arg  => undef,
		  default   => sub{ Xray::FeffPathWrapper->new() },
		  #lazy      => 1,
		  #builder   => '_build_object',
		 );

has 'errorcode'    => (is => 'rw', isa => 'Int',  default => 0,);
has 'errormessage' => (is => 'rw', isa => 'Str',  default => q{},);


has 'phpad'   => (is => 'rw', isa => 'Str',  default => 'phase.pad',    trigger => sub{pushback(@_, 'phpad'  )});

has 'Index'   => (is => 'rw', isa => 'Int',  default => 9999, trigger => sub{pushback(@_, 'Index'  )},);
has 'nleg'    => (is => 'rw', isa => 'Int',  default => 0,    trigger => sub{pushback(@_, 'nleg'   )},);
has 'degen'   => (is => 'rw', isa => 'Num',  default => 1.0,  trigger => sub{pushback(@_, 'degen'  )},);
has 'iorder'  => (is => 'rw', isa => 'Int',  default => 2,    trigger => sub{pushback(@_, 'iorder' )},);

has 'nnnn'    => (is => 'rw', isa => 'Bool', default => 0,    trigger => sub{pushback(@_, 'nnnn'   )},);
has 'json'    => (is => 'rw', isa => 'Bool', default => 0,    trigger => sub{pushback(@_, 'json'   )},);
has 'verbose' => (is => 'rw', isa => 'Bool', default => 0,    trigger => sub{pushback(@_, 'verbose')},);

has 'ipol'    => (is => 'rw', isa => 'Int',  default => 0,    trigger => sub{pushback(@_, 'ipol'   )},);
has 'elpty'   => (is => 'rw', isa => 'Num',  default => 0.0,  trigger => sub{pushback(@_, 'elpty'  )},);

has 'edge'    => (is => 'rw', isa => 'Num',  default => 0.0);
has 'gam_ch'  => (is => 'rw', isa => 'Num',  default => 0.0);
has 'kf'      => (is => 'rw', isa => 'Num',  default => 0.0);
has 'mu'      => (is => 'rw', isa => 'Num',  default => 0.0);
has 'rnorman' => (is => 'rw', isa => 'Num',  default => 0.0);
has 'rs_int'  => (is => 'rw', isa => 'Num',  default => 0.0);
has 'vint'    => (is => 'rw', isa => 'Num',  default => 0.0);

has 'exch'    => (is => 'rw', isa => 'Str',  default => '');
has 'version' => (is => 'rw', isa => 'Str',  default => '');

has 'absorber' => (is => 'rw', isa => 'ArrayRef', default => sub{[0,0,0]});

## arrays
has 'evec'    => (traits  => ['Array'],
                  is      => 'rw',
                  isa     => 'ArrayRef[Num]',
                  default => sub { [0,0,0] },
                  trigger => \&evec_set,
                 );
has 'xivec'   => (traits  => ['Array'],
                  is      => 'rw',
                  isa     => 'ArrayRef[Num]',
                  default => sub { [0,0,0] },
                  trigger => \&xivec_set);

has 'ipot'     => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Int]', default => sub { [] },
                   handles => {clear_ipot => 'clear', push_ipot  => 'push', });
has 'rat'      => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[ArrayRef[Num]]', default => sub { [] },
                   handles => {clear_rat => 'clear', push_rat  => 'push', });
has 'iz'       => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Int]', default => sub { [] },
                   handles => {clear_iz => 'clear', push_iz  => 'push', });

has 'reff'     => (is => 'rw', isa => 'Num',  default => 0.0);
has 'ri'       => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
                   handles => {clear_ri => 'clear', push_ri  => 'push', });
has 'beta'     => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
                   handles => {clear_beta => 'clear', push_beta  => 'push', });
has 'eta'      => (traits  => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub { [] },
                   handles => {clear_eta => 'clear', push_eta  => 'push', });





has 'ne'       => (is => 'rw', isa => 'Int',  default => 0);
has 'k'        => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'real_phc' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'mag_feff' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'pha_feff' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'red_fact' => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'lam'      => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'rep'      => (is => 'rw', isa => 'ArrayRef', default => sub{[]});

  # foreach my $att ($self->meta->get_all_attributes) {
  #   print $att->name, $/;
  # };

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


sub pushback {
  my ($self, $new, $old, $which) = @_;
  return if (any {$_ eq $which} qw(nleg ne edge gam_ch kf mu rnorman version exch rs_int vint));
  my $method = '_set_' . lc($which);
  if ($self->meta->get_attribute($which)->type_constraint eq 'Num') {
    $self->wrapper->$method(1.0*$new);
  } elsif ($self->meta->get_attribute($which)->type_constraint eq 'Bool') {
    my $val = ($new) ? 1 : 0;
    $self->wrapper->$method($val);
    $self->phpad($self->phpad);	# what the hell?  why does this get unset when setting a boolean????
  } else {
    $self->wrapper->$method($new);
  };
  $self->ipol(1) if (($which eq 'elpty') and $self->wrapper->_elpty);
};


sub clear {
  my ($self) = @_;
  $self->wrapper->_clear_path;
  foreach my $att ($self->meta->get_all_attributes) {
    next if ($self->meta->get_attribute($att->name)->type_constraint =~ m{FeffPath});
    my $attr = $att->name;
    my $method = '_'.lc($att->name);
    if ($attr =~ m{(?:(?:e|xi)vec|absorber)}) {
      $self->$attr([0,0,0])
    } elsif ($self->meta->get_attribute($att->name)->type_constraint =~ m{ArrayRef}) {
      $self->$attr([])
    } elsif ($self->meta->get_attribute($att->name)->type_constraint =~ m{Bool}) {
      my $val = $self->wrapper->$method || 0;
      $self->$attr(0);
    } else {
      my $val = $self->wrapper->$method || 0;
      $self->$attr($val);
    };
  };
  return $self;
};

sub evec_set {
  1;
};
sub xivec_set {
  1;
};

sub atom {
  my ($self, $x, $y, $z, $ip) = @_;
  my $message;
  my $xx = $x - $self->absorber->[0];
  my $yy = $y - $self->absorber->[1];
  my $zz = $z - $self->absorber->[2];
  my $err = $self->wrapper->_add_scatterer($xx, $yy, $zz, $ip);
  $self->errorcode($err);
  my $em = $self->wrapper->_errormessage;
  $em =~ s{add_scatterer}{atom method};
  $self->errormessage($em);
  $self->nleg($self->wrapper->_nleg);
  $self->push_ipot($ip);
  $self->push_rat([$x, $y, $z]);
  return $err;
};

sub path {
  my ($self) = @_;
  my $message;
  $self->push_ipot(0);
  $self->push_rat($self->absorber);
  my $err = $self->wrapper->_make_path;
  $self->errorcode($err);
  my $em = $self->wrapper->_errormessage;
  $em =~ s{make_path}{path method};
  $self->errormessage($em);
  if (not $err) {
    $self->ne($self->wrapper->_ne);
    foreach my $a (qw(k real_phc mag_feff pha_feff red_fact lam rep)) {
      my $method = '_' . $a . '_array';
      my @array = $self->wrapper->$method;
      $self->$a(\@array);
    };
    foreach my $a (qw(ri beta eta iz)) {
      my $method = '_' . $a . '_array';
      my @array = $self->wrapper->$method;
      $self->$a(\@array);
    };
  };
};


no Moose;
__PACKAGE__->meta->make_immutable;
1;



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

=head1 BUGS AND LIMITATIONS

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
Last update: 4 November, 2014

=cut
