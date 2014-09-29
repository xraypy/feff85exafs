package Xray::FeffPath;

BEGIN {
  $ENV{FEFFSRC} = '/home/bruce/git/feff85exafs/src/';
}


use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
# our @EXPORT_OK = (  );
# our @EXPORT = qw(  );

our $VERSION = '1.00';

use Inline C => 'DATA',
           LIBS => '-L'.$ENV{FEFFSRC}.'GENFMT -L'.$ENV{FEFFSRC}.'COMMON -L'.$ENV{FEFFSRC}.'MATH -L'.$ENV{FEFFSRC}.'PAR -lfeffpath -lonepath -lgenfmt -lfeffcom -lfeffmath -lfeffpar -lgfortran -lm',
           NAME => 'Xray::FeffPath';

1;
__DATA__

=pod


=head1 NAME

Xray::FeffPath - Inline::C interface to libfeffpath

=head1 VERSION

2.00

=head1 SYNOPSIS

This provides a simple object oriented interface to the C<feffpath>
wrapper around Feff's C<onepath> subroutine.  C<onepath> combines the
F-matrix calculation in F<GENFMT/genfmt.f> with the construction of the
common presentation of F-effective from F<FF2X/feffdt.f>.

Given the geometry of a scattering path in the form of Cartesian
coordinates of atoms making up the path, along with some additional
information such as the degeneracy, compute F-effective for that
scattering path.

The following computes the columns of the F<feff0004.dat> file for
copper metal:

  #!/usr/bin/perl
  use strict;
  use warnings;
  use Xray::FeffPath;

  my $errcode;
  my $path = Xray::FeffPath->new($errcode);
  $path->atom(0, 0, -3.61, 1);
  $path->atom(-1.805, 0, -1.805, 1);
  $path->path;
  my @kgrid  = $path->kgrid;
  my @caps   = $path->real_phc;
  my @amff   = $path->mag_feff;
  my @phff   = $path->pha_feff;
  my @redfac = $path->red_fact;
  my @lambda = $path->lam;
  my @rep    = $path->rep;

This closely follows the example at
L<https://metacpan.org/module/Inline::C-Cookbook#Object-Oriented-Inline>,
mapping all the functions demonstrated in F<xdi_reader.c>.

See L<Xray::XDI> for a Moose-ified wrapper around this.

All methods share a name with L<Xray::XDI>, except that methods of
this object are all preceded with an underscore.

=head1 METHODS

=head2 action methods

=over 4

=item C<new>

Create the FeffPath object

   my $errcode;
   my $path = Xray::FeffPath->new($errcode);

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

=head2 getter methods

Each of the members of the FEFFPATH struct is wrapped in a getter
method of the same name.  Each of the following exists:

    my $path = Xray::FeffPath->new($errcode);
    $index   = $path->Index;    # path index
    $deg     = $path->deg;      # path degeneracy
    $nleg    = $path->nleg;     # number of legs
    $iorder  = $path->iorder;   # order of MS approximation
    $nnnn    = $path->nnnn;     # flag for writing feffNNNN.dat
    $json    = $path->json;     # flag for writing feffNNNN.json
    $verbose = $path->verbose;  # flag for writing screen messages
    $ipol    = $path->ipol;     # flag for doing polarization calculation
    @evec    = $path->evec;     # polarization vector
    $elpty   = $path->elpty;    # ellipticity
    @xivec   = $path->xivec;    # poynting vector
    @ri      = $path->ri;       # leg lengths
    @beta    = $path->beta;     # beta angles
    @eta     = $path->eta;      # eta angles
    $reff    = $path->reff;     # half path length
    $ne      = $path->ne;       # length of k grid
    @kgrid   = $path->kgrid;    # column 1 from feffNNNN.dat
    @caps    = $path->real_phc; # column 2 from feffNNNN.dat
    @amff    = $path->mag_feff; # column 3 from feffNNNN.dat
    @phff    = $path->pha_feff; # column 4 from feffNNNN.dat
    @redfac  = $path->red_fact; # column 5 from feffNNNN.dat
    @lam     = $path->lam;      # column 6 from feffNNNN.dat
    @rep     = $path->rep;      # column 7 from feffNNNN.dat

Note that the getter method for the path index is capitalized.  Also
note that getter methods for array members of the FEFFPATH struct
return arrays, not array references.  C<@evec> and C<@xivec> are 3
elements long.  C<@ri>, C<@beta>, and C<@eta> are C<$nleg> elements
long.  The arrays for the 7 columns of F<feffNNNN.dat> are C<$ne>
elements long.

=head2 setter methods

For members of the FEFFPATH struct that make sense to set from a
script, setter methods exist:

    my $path = Xray::FeffPath->new($errcode);
    $path->set_index($val);     # path index
    $path->set_deg($val);       # path degeneracy
    $path->set_iorder($val);    # order of MS approximation
    $path->set_nnnn($val);      # flag for writing feffNNNN.dat
    $path->set_json($val);      # flag for writing feffNNNN.json
    $path->set_verbose($val);   # flag for writing screen messages
    $path->set_ipol($val);      # flag for doing polarization calculation
    $path->set_evec($x,$y,$z);  # polarization vector
    $path->set_elpty($val);     # flag for doing polarization calculation
    $path->set_xivec($x,$y,$z); # poynting vector

=head1 DEPENDENCIES

=over 4

=item *

L<Inline>

=item *

L<Inline::C>

=back

=head1 BUGS AND LIMITATIONS

The C section of this file is mostly cargo-cult programming taken from
the Object Oriented Inline example at
L<https://metacpan.org/pod/Inline::C::Cookbook>.  The definition of
Newx should have this comment block explaining it:

    Allocate memory with Newx if it's available - if it's an older
    perl that doesn't have Newx then we resort to using New.

Regarding INT2PTR, see
Lhttp://www.mail-archive.com/inline%40perl.org/msg02689.html>.  The
use of PTR2IV in new was the result of trial and error, rather than
deep understanding.

=over 4

=item *

Getter method(s) for rat and ipot

=item *

F<feffpath.h> needs to be installed someplace findable

=item *

The parts of feff need to be installed somewhare that the C<LIBS>
argument to Inline can find them.

=item *

Error handling is nonexistent

=back

Please report problems to Bruce Ravel (bravel AT bnl DOT gov)

Patches are welcome.

=head1 AUTHOR

Bruce Ravel (bravel AT bnl DOT gov)

L<http://cars9.uchicago.edu/~ravel/software/>

=head1 LICENCE AND COPYRIGHT

To the extent possible, the authors have waived all rights granted by
copyright law and related laws for the code and documentation that
make up the Perl Interface to the XAS Data Interchange Format.
While information about Authorship may be retained in some files for
historical reasons, this work is hereby placed in the Public Domain.
This work is published from: United States.

Note that the feffpath library itself is NOT public domain.

Author: Bruce Ravel (bravel AT bnl DOT gov).
Last update: 25 September, 2014

=cut


__C__

#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif

#include "feffpath.h"

SV* new(char* class) {
  FEFFPATH * path;
  int        ret;
  SV *       obj;
  SV *       obj_ref;

  Newx(path, 1, FEFFPATH);

  ret = create_path(path);
  /*sv_setiv(errcode, ret);*/

  obj = newSViv((IV)PTR2IV(path));
  obj_ref = newRV_noinc(obj);

  sv_bless(obj_ref, gv_stashpv(class, GV_ADD));
  SvREADONLY_on(obj);
  return obj_ref;
}

int atom(SV* obj, double xx, double yy, double zz, int ip) {
  int nleg;
  nleg = add_scatterer((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))), xx, yy, zz, ip);
  return nleg;
}

int path(SV* obj) {
  int ret;
  ret = make_path((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))));
  return ret;
};

void clear(SV* obj) {
  int ret;
  clear_path((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))));
};


bool nnnn(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nnnn;
}

bool json(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->json;
}

bool verbose(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->verbose;
}

int Index(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->index;
}

int nleg(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg;
}

int ipol(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ipol;
}

int iorder(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->iorder;
}

int ne(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne;
}

double deg(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->deg;
}

double elpty(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->elpty;
}

double reff(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->reff;
}

/* polarization and ellipticity vectors */

void evec(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < 3; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->evec[i])));
  }
  Inline_Stack_Done;
}

void xivec(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < 3; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->xivec[i])));
  }
  Inline_Stack_Done;
}

/* path geometry */

void ri(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ri[i])));
  }
  Inline_Stack_Done;
}

void beta(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->beta[i])));
  }
  Inline_Stack_Done;
}

void eta(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->eta[i])));
  }
  Inline_Stack_Done;
}

/* columns from feffNNNN.dat */

void kgrid(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->k[i])));
  }
  Inline_Stack_Done;
}

void real_phc(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->real_phc[i])));
  }
  Inline_Stack_Done;
}

void mag_feff(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->mag_feff[i])));
  }
  Inline_Stack_Done;
}

void pha_feff(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->pha_feff[i])));
  }
  Inline_Stack_Done;
}

void red_fact(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->red_fact[i])));
  }
  Inline_Stack_Done;
}

void lam(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->lam[i])));
  }
  Inline_Stack_Done;
}

void rep(SV* obj) {
  long i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->rep[i])));
  }
  Inline_Stack_Done;
}



/* set boolean valued input to the FEFFPATH struct */

void set_nnnn(SV* obj, int val) {
  bool bval;
  bval = (val == 0) ? false : true;
  ((FEFFPATH*)SvIV(SvRV(obj)))->nnnn = bval;
}

void set_json(SV* obj, int val) {
  bool bval;
  bval = (val == 0) ? false : true;
  ((FEFFPATH*)SvIV(SvRV(obj)))->json = bval;
}

void set_verbose(SV* obj, int val) {
  bool bval;
  bval = (val == 0) ? false : true;
  ((FEFFPATH*)SvIV(SvRV(obj)))->verbose = bval;
}

/* set path geometry scalars in the FEFFPATH struct */

void set_index(SV* obj, int val) {
  ((FEFFPATH*)SvIV(SvRV(obj)))->index = val;
}

void set_iorder(SV* obj, int val) {
  ((FEFFPATH*)SvIV(SvRV(obj)))->iorder = val;
}

void set_deg(SV* obj, double val) {
  ((FEFFPATH*)SvIV(SvRV(obj)))->deg = val;
}

/* set polarization and ellipticity members */

void set_ipol(SV* obj, int val) {
  ((FEFFPATH*)SvIV(SvRV(obj)))->ipol = val;
}

void set_elpty(SV* obj, double val) {
  ((FEFFPATH*)SvIV(SvRV(obj)))->elpty = val;
}

void set_evec(SV* obj, double xx, double yy, double zz) {
  ((FEFFPATH*)SvIV(SvRV(obj)))->evec[0] = xx;
  ((FEFFPATH*)SvIV(SvRV(obj)))->evec[1] = yy;
  ((FEFFPATH*)SvIV(SvRV(obj)))->evec[2] = zz;
}

void set_xivec(SV* obj, double xx, double yy, double zz) {
  ((FEFFPATH*)SvIV(SvRV(obj)))->xivec[0] = xx;
  ((FEFFPATH*)SvIV(SvRV(obj)))->xivec[1] = yy;
  ((FEFFPATH*)SvIV(SvRV(obj)))->xivec[2] = zz;
}


void DESTROY(SV* obj) {
  FEFFPATH* path = (FEFFPATH*)SvIV(SvRV(obj));
  cleanup(path);
}

