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

use Term::ANSIColor qw(:constants);
sub trace {
  my ($self) = @_;
  my $max_depth = 30;
  my $i = 0;
  my $base = substr($INC{'Demeter.pm'}, 0, -10);
  my ($green, $red, $yellow, $end) = (BOLD.GREEN, BOLD.RED, BOLD.YELLOW, RESET);
  local $|=1;
  print($/.BOLD."--- Begin stack trace ---$end\n");
  while ( (my @call_details = (caller($i++))) && ($i<$max_depth) ) {
    (my $from = $call_details[1]) =~ s{$base}{};
    my $line  = $call_details[2];
    my $color = RESET.YELLOW;
    (my $func = $call_details[3]) =~ s{(?<=::)(\w+)\z}{$color$1};
    print("$green$from$end line $red$line$end in function $yellow$func$end\n");
  }
  print(BOLD."--- End stack trace ---$end\n");
  return $self;
};



no Moose;
# no need to fiddle with inline_constructor here
__PACKAGE__->meta->make_immutable;
