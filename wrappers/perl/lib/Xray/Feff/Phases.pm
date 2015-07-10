package Xray::Feff::Phases;

use Moose;
#use MooseX::Storage;
#with Storage('format' => 'JSON', 'io' => 'File');
use MooseX::NonMoose;
use MooseX::Aliases;
extends 'Xray::Feff::PhasesWrapper';


use JSON;
use List::MoreUtils qw(any);

our $VERSION = '1.00'; # Inline::MakeMake uses /^\d.\d\d$/ as the
                       # pattern for the version number -- note the
                       # two digits to the right of the dot

has 'wrapper' => (
		  is        => 'ro',
		  #traits => [qw(DoNotSerialize)],
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

## i/o files
has 'phpad'    => (is => 'rw', isa => 'Str',  default => 'phase.pad',     trigger => sub{pushback(@_, 'phpad'    )});
has 'jsonfile' => (is => 'rw', isa => 'Str',  default => 'libpotph.json', trigger => sub{pushback(@_, 'jsonfile' )});
has 'useperljson' => (is => 'rw', isa => 'Bool',  default => 1);
has 'verbose'  => (is => 'rw', isa => 'Bool',  default => 0,              trigger => sub{pushback(@_, 'verbose'  )});

## scalars
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
has 'gamach'   => (is => 'rw', isa => 'Num',  default => 0.0,             trigger => sub{pushback(@_, 'gamach'   )});
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

## three vectors
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
#has 'ptz'      => (is => 'rw', isa => 'ArrayRef[ArrayRef]', default => sub{[]});


## potentials information
has 'iz'       => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "iz")});
has 'potlbl'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "potlbl")});
has 'lmaxsc'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "lmaxsc")});
has 'lmaxph'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "lmaxph")});
has 'xnatph'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "xnatph")});
has 'spinph'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "spinph")});
has 'folp'     => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "folp")});
has 'xion'     => (is => 'rw', isa => 'ArrayRef', default => sub{[]},  trigger => sub{set_pot_array(@_, "xion")});

## atoms list
has 'iphat'    => (is => 'rw', isa => 'ArrayRef', default => sub{[]});
has 'rat'      => (is => 'rw', isa => 'ArrayRef[ArrayRef]', default => sub{[]});

has 'titles'   => (is => 'rw', isa => 'ArrayRef', default => sub{[]}, trigger => \&set_titles_array);



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

sub clear {
  my ($self) = @_;
  $self->wrapper->_clear_phases;
  return $self;
};

sub _dump {
  my ($self) = @_;
  $self->wrapper->_dump;
  return $self;
};

sub _json_c {
  my ($self) = @_;
  my $ret = $self->wrapper->_read_libpotph_json;
  $self->_fetch;
  return $ret;
};
sub _json_pp {
  my ($self) = @_;
  my $json_text = $self->_slurp($self->jsonfile);
  my $rhash = decode_json($json_text);

  foreach my $meth (qw(ntitle nat nph ihole rscf lscf nscmt ca nmix ecv icoul ipol elpty ispin angks gamach
		      ixc vr0 vi0 ixc0 iafolp rgrd iunf inters totvol jumprm nohole
		      iz potlbl lmaxsc lmaxph xnatph spinph folp xion titles)) {
    my $att = ($meth eq 'nat')   ? 'natt'
            : ($meth eq 'rscf')  ? 'rfms1'
            : ($meth eq 'lscf')  ? 'lfms1'
            : ($meth eq 'ca')    ? 'ca1'
            : ($meth eq 'vr0')   ? 'vro'
            : ($meth eq 'vi0')   ? 'vio'
            : ($meth eq 'iphat') ? 'iphatx'
            : $meth;
    $self->$meth($rhash->{$att});
  };
  ## truncate lists to correct size
  foreach my $meth (qw(iz potlbl lmaxsc lmaxph xnatph spinph folp xion)) {
    $#{$self->$meth} = $self->nph;
  };
  $#{$self->titles} = $self->ntitle;
  my @lables = ();
  foreach my $pl (@{$self->potlbl}) {
    push @lables, sprintf("%-6s", $pl);
  };
  $self->potlbl(\@lables);
  my @titles = ();
  foreach my $ti (@{$self->titles}) {
    push @titles, sprintf("%-79s\0", $ti);
  };
  $self->titles(\@titles);

  my @rat = ();
  my $x   = $rhash->{x};
  my $y   = $rhash->{y};
  my $z   = $rhash->{z};
  foreach my $i (0 .. $rhash->{natt}-1) {
    my $this = [$x->[$i], $y->[$i], $z->[$i]];
    push @rat, $this;
  };
  $self->rat(\@rat);
  $self->wrapper->_set_rat_array(map {@$_} @rat); # flatten argument, then unflatten for the C struct

  my $ip = $rhash->{iphatx};
  $self->iphat($ip);
  $self->wrapper->_set_iphat_array(@$ip);
};


sub phases {
  my ($self) = @_;
  my $ret = $self->wrapper->_make_phases;
  ## fetch polarization tensor
  $self->errormessage($self->wrapper->_errormessage);
  $self->errorcode($self->wrapper->_errorcode);
  if ($ret) {
    printf("%s\n\t(error code %d)\n", $self->errormessage, $self->errorcode);
  };
  return $ret;
};

sub _fetch {
  my ($self) = @_;
  foreach my $att (qw(ntitle nat nph ihole rscf lscf nscmt ca nmix ecv icoul ipol elpty ispin angks gamach
		      ixc vr0 vi0 ixc0 iafolp rgrd iunf inters totvol jumprm nohole)) {
    my $method = '_'.$att;
    $self->$att($self->wrapper->$method);
  };
  foreach my $att (qw(evec xivec spvec)) {
    my $method = '_'.$att;
    $self->$att([$self->wrapper->$method]);
  };
  foreach my $att (qw(iz potlbl lmaxsc lmaxph xnatph spinph folp xion titles iphat)) { ## rat
    my $method = '_'.$att.'_array';
    $self->$att([$self->wrapper->$method]);
  };
  my @biglist = $self->wrapper->_rat_array;
  my @rat = ();
  while (@biglist) {		# rat[nat][3] is returned as a flattened list (x1,y1,z1,x2,y2,z2,...xN,yN,zN)
    my $this = [shift(@biglist), shift(@biglist), shift(@biglist)];
    push @rat, $this;
  };
  $self->rat(\@rat);
  return $self;
};


## ------------------------------------------------------------
## trigger methods

sub pushback {
  my ($self, $new, $old, $which) = @_;
  #print "setting $which\n";
  #return if (any {$_ eq $which} qw(gamach));
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
  if ($which eq 'jsonfile') {
    ($self->useperljson) ? $self->_json_pp : $self->_json_c;
  };
};

sub evec_set {
  my ($self, $vec) = @_;
  $self->wrapper->_set_evec(1.0*$vec->[0], 1.0*$vec->[1], 1.0*$vec->[2]);
  $self->ipol(1) if ($vec->[0] or $vec->[1] or $vec->[2]);
  return $self;
};
sub xivec_set {
  my ($self, $vec) = @_;
  $self->wrapper->_set_xivec(1.0*$vec->[0], 1.0*$vec->[1], 1.0*$vec->[2]);
  $self->ipol(1) if ($vec->[0] or $vec->[1] or $vec->[2]);
  return $self;
};
sub spvec_set {
  my ($self, $vec) = @_;
  $self->wrapper->_set_spvec(1.0*$vec->[0], 1.0*$vec->[1], 1.0*$vec->[2]);
  $self->ispin(1) if ($vec->[0] or $vec->[1] or $vec->[2]);
  return $self;
};

sub set_pot_array {
  my ($self, $new, $old, $which) = @_;
  #print "setting $which\n";
  my @list;
  if ($which =~ m{iz|iphat|lmaxsc|lmaxph}) {
    @list = map {int($_)} @$new; # constrain to be integer
    #print join("|", @list), $/ if $which eq 'iz';
  } elsif ($which =~ m{xnatph|folp|xion|spinph}) {
    @list = map {$_*1.0} @$new;	# constrain to be float
  } elsif ($which =~ m{potlbl}) {
    foreach my $s (@$new) { # constrain to be 6 characters or less
      push @list, (length($s) > 6) ? substr($s,0,6) : $s;
    };
  };
  my $method = '_set_' . lc($which) . '_array';
  $self->wrapper->$method(@list);
  return $self;
};

sub set_titles_array {
  my ($self, $new, $old) = @_;
  my @list;
  foreach my $s (@$new) { # constrain to be 6 characters or less
    push @list, (length($s) > 6) ? substr($s,0,6) : $s;
  };
  $self->wrapper->_set_titles_array(@list);
  return $self;
};

sub set_rat_array {
  my ($self, $new, $old) = @_;
  my @list = map {@$_} @{$new}; # flatten rat, which is a list of lists
  $self->wrapper->_set_rat_array(@list);
  return $self;
};

sub _slurp {
  my ($class, $file) = @_;
  local $/;
  return q{} if (not -e $file);
  return q{} if (not -r $file);
  open(my $FH, $file);
  my $text = <$FH>;
  close $FH;
  return $text;
};

sub serialize {
  my ($self, $fname) = @_;
  $self->store($fname);
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
  my $cu = Xray::Feff::Phases->new;
     ##--- gather information about the calculation and cluster
  $cu->jsonfile('../fortran/libpotph.json');
     ##--- set the output file
  $cu->phpad('cu_phases.pad');
     ##--- run feff's potentials and phases calculator
  $cu->phases;

=head1 INSTALLATION

See L<Xray::Feff>.

=head1 METHODS

=over 4

=item C<new>

Create the FeffPhases object

   my $phases = Xray::Feff::Phases->new();

=item C<clear>

Reinitialize the phases object

   $path->clear;

=item C<phases>

Compute the phases and write the output file, F<phase.pad> or as
specified by the C<phpad> attribute.

   $path->phases;

=back

To destroy a FeffPath object, simply do

   undef $path;

Deallocation of memory in the C stuct will be handled by the wrapper.

=head1 ATTRIBUTES

Each of the members of the FEFFPATH struct is wrapped in a getter
method of the same name.  Each of the following exists:

=head2 Scalar valued attributes

=over 4

=item C<phpad> (character, default = phase.pad)

The path for writing the F<phase.pad> file.

=item C<jsonfile> (character, default = libpotph.json)

The path to the F<libpth.json> file.  The JSON file is read
immediately upon setting this attribute.

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

=iten C<verbose> (boolean)

Controls whether Feff writes its screen messages.  False means to run
Feff silently.

=item C<ntitle> (int)

number of header lines (TITLE)

=item C<titles> (\*\*char)

(nheadx) array of header string (TITLE)

=item C<nat> (int)

number of atoms in cluster (ATOMS)

=item C<rat> (\*\*double)

(3,natx) cartesian coordinates of atoms in cluster (ATOMS)

=item C<iphat> (\*int)

(natx) unique potential indeces of atoms in cluster (ATOMS)

=item C<nph> (int)

number of unique potentials (POTENTIALS)

=item C<iz> (\*int)

(0:nphx) Z numbers of unique potentials (POTENTIALS)

=item C<potlbl> (\*\*char)

(0:nphx) labels of unique potentials (POTENTIALS)

=item C<lmaxsc> (\*int)

(0:nphx) l max for SCF for each potential (POTENTIALS)

=item C<lmaxph> (\*int)

(0:nphx) l max for FMS for each potential (POTENTIALS)

=item C<xnatph> (\*double)

(0:nphx) stoichiometry of each potential (POTENTIALS)

=item C<spinph> (\*double)

(0:nphx) spin on each unique potential (POTENTIALS)

=item C<ihole> (int)

edge index, 1=K, 4=L3, etc (HOLE/EDGE)

=item C<rscf> (double)

cluster radius for self-consistent calculation (SCF)

=item C<lscf> (int)

0=solid, 1=molecule (SCF)

=item C<nscmt> (int)

max number of self-consistency iterations (SCF)

=item C<ca> (double)

self-consistency convergence accelerator (SCF)

=item C<nmix> (int)

number of mixing iterations before Broyden (SCF)

=item C<ecv> (double)

core/valence separation energy (SCF)

=item C<icoul> (int)

obsolete param. for handling Coulomb potential (SCF)

=item C<ipol> (int)

1=do polarization calculation (POLARIZATION)

=item C<evec> (\*double)

(3) polarization array (POLARIZATION)

=item C<elpty> (double)

eccentricity of elliptical light (ELLIPTICITY)

=item C<xivec> (\*double)

(3) ellipticity array (ELLIPTICITY)

=item C<ispin> (int)

+/-2 = do spin calculation (SPIN)

=item C<spvec> (\*double)

(3) spin array (SPIN)

=item C<angks> (double)

angle between spin and incidient beam (SPIN)

=item C<ptz> (double complex)

(-1:1,-1:1) polarization tensor (return)

=item C<gamach> (double)

tabulated core-hole lifetime (return)

=item C<ixc> (int)

exchange index (EXCHANGE)

=item C<vr0> (double)

Fermi level offset (EXCHANGE)

=item C<vi0> (double)

constant broadening (EXCHANGE)

=item C<ixc0> (int)

exchange index for background function (EXCHANGE)

=item C<iafolp> (int)

1=do automated overlapping (FOLP & AFOLP)

=item C<folp> (\*double)

(0:nphx) overlapping fractions (FOLP & AFOLP)

=item C<xion> (\*double)

(0:nphx) potential ionizations (ION)

=item C<rgrd> (double)

radial grid used for the potentials/phases (RGRID)

=item C<iunf> (int)

1=unfreeze f electrons (UNFREEZEF)

=item C<inters> (int)

(INTERSTITIAL)

=item C<totvol> (double)

(INTERSTITIAL)

=item C<jumprm> (int)

1=remove potential jumps at muffin tin radii (JUMPRM)

=item C<nohole> (int)

1=compute without core-hole (NOHOLE)

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

=head2 3-vector array-reference

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

=back

=head2 Array-references decribing the unique potentials

=over 4

=item C<iz>

A reference to an array containing the Z numbers of each unique
potential.

=item C<potlbl>

Label for each potential index (up to 6 characters)

=item C<lmaxsc>

Maximum angular momentum for use in self consistency calculation

=item C<lmaxph>

Maximum angular momentum for use in XANES calculation

=item C<xnatph>

Stoichiometry of each potential in the cluster

=item C<spinph>

Spin on each potential

=item C<folp>

Overlap fraction for each potential

=item C<xion>

Ionization on each potential.

=back

=head2 Array-references decribing the unique potentials

=over 4

=item C<iphat>

=item C<rat>

=back

=head1 Error codes

The error codes returned by the C<phases> methods are the sums of the
codes actually found by the method call, thus the returned error codes
are meant to be interpreted bitwise.

=over 4

=item I<1>

Too many unique potentials

=item I<2>

Too many atoms

=item I<4>

Edge index ust be between 1 and 9, i.e. K to M5

=item I<8>

A potential is defined with an invalid Z number

=item I<16>

A potential is defined with an invalid angular momentum

=item I<32>

A potential is defined with a negative stoichiometry

=item I<64>

Overlap fraction must be between 0.7 and 1.5

=item I<128>

invalid ionization in potentials list (NOT USED)

=item I<256>

bad rscf (NOT USED)

=item I<512>

Convergence accelerator (ca) should be around 0.2, maybe a bit smaller

=item I<1024>

Core/valence separation energy (ecv) must be negative

=item I<2048>

Exchange index (ixc) be 0, 1, 2, 3, or 5

=item I<4096>

Radial grid (rgrid) should be around 0.05, maybe a bit smaller, not negative

=item I<8192>

Invalid potential index used or additional absorber in cluster

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

Need better methods for moving C<rat> and C<iphat> into and out of the
C struct.  I think I need to know a bit more about the eventual UI
intended to replace the F<feff.inp> file....

=item *

Polarization tensor is not handled at all.

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
