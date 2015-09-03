package Xray::Feff::PathWrapper;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
# our @EXPORT_OK = (  );
# our @EXPORT = qw(  );

our $VERSION = '1.00';

use Inline C => 'DATA',
           LIBS => '-lfeffpath',
           VERSION => '1.00',
           NAME => 'Xray::Feff::PathWrapper';


## for Windows, will need
##   LIBS => '-LC:\\Program Files\\larch\\dlls\\win32 -LC:\\Program Files\\larch\\dlls\\win64 -lfeffpath'
##   INC  => '-IC:\\Program Files\\larch\\include'
##
## also will need Strawberry with gcc 4.9.  Yipes!

1;
__DATA__



=head1 NAME

Xray::FeffPath::Wrapper - Inline::C interface to libfeffpath


=head1 VERSION


=head1 SYNOPSIS


=head1 DEPENDENCIES

=over 4

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

=head1 AUTHOR

Bruce Ravel (bravel AT bnl DOT gov)

L<http://cars9.uchicago.edu/~ravel/software/>

=head1 LICENCE AND COPYRIGHT

To the extent possible, the authors have waived all rights granted by
copyright law and related laws for the code and documentation that
make up the Perl Interface to the FeffPath library.  While information
about Authorship may be retained in some files for historical reasons,
this work is hereby placed in the Public Domain.  This work is
published from: United States.

Note that the feffpath and onepath libraries themselves are NOT public
domain, nor is the Fortran source code for Feff that it relies upon.

Author: Bruce Ravel (bravel AT bnl DOT gov).
Last update: 13 February, 2015

=cut

__C__

#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif

#include <string.h>
#include "feffpath.h"

SV* new(char* class) {
  FEFFPATH * feffpathwrapper;
  SV *      obj;
  SV *      obj_ref;

  Newx(feffpathwrapper, 1, FEFFPATH);

  obj = newSViv((IV)PTR2IV(feffpathwrapper));
  obj_ref = newRV_noinc(obj);

  sv_bless(obj_ref, gv_stashpv(class, GV_ADD));
  SvREADONLY_on(obj);
  return obj_ref;
}

int    _nex   (SV* obj) { return nex;    }
int    _nphx  (SV* obj) { return nphx;   }
int    _npatx (SV* obj) { return npatx;  }
int    _legtot(SV* obj) { return legtot; }
double _bohr  (SV* obj) { return bohr;   }
double _ryd   (SV* obj) { return ryd;    }
double _hart  (SV* obj) { return hart;   }

int _create_path(SV* obj) {
  int ret;
  ret = create_path((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))));
  return ret;
}
void _clear_path(SV* obj) {
  clear_path((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))));
}
void _cleanup(SV* obj) {
  cleanup((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))));
}
int _add_scatterer(SV* obj, double xx, double yy, double zz, int ip) {
  int ret;
  ret = add_scatterer((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))), xx, yy, zz, ip);
  return ret;
}
int _make_path(SV* obj) {
  int ret;
  ret = make_path((INT2PTR(FEFFPATH*, SvIV(SvRV(obj)))));
  return ret;
}


int _errorcode(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->errorcode;
}
char* _errormessage(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->errormessage;
}



/* scalar valued attributes */

char* _phpad(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->phpad;
}
void _set_phpad(SV* obj, char* c) {
       strcpy((INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->phpad, c);
}

/* ----- path index */
int _index(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->index;
}
void _set_index(SV* obj, int i) {
       (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->index = i;
}

/* ----- number of legs */
int _nleg(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg;
}

/* ----- degeneracy */
double _degen(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->degen;
}
void _set_degen(SV* obj, double d) {
       (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->degen = d;
}

/* ----- iorder */
int _iorder(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->iorder;
}
void _set_iorder(SV* obj, int i) {
       (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->iorder = i;
}

/* ----- flag to write feffNNNN.dat */
bool _nnnn(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nnnn;
}
void _set_nnnn(SV* obj, bool i) {
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nnnn = i;
}

/* ----- flag to write feffNNNN.xdi */
bool _xdi(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->xdi;
}
void _set_xdi(SV* obj, bool i) {
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->xdi = i;
}

/* ----- flag to be verbose */
bool _verbose(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->verbose;
}
void _set_verbose(SV* obj, bool b) {
       (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->verbose = b;
}

/* ----- flag for polarization & ellipticity calculation */
bool _ipol(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ipol;
}
void _set_ipol(SV* obj, bool b) {
       (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ipol = b;
}

/* ----- ellipticity value */
double _elpty(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->elpty;
}
void _set_elpty(SV* obj, double d) {
       (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->elpty = d;
}

/* ----- 3-vector valued attributes */

void _evec(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < 3; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->evec[i])));
  }
  Inline_Stack_Done;
}
void _set_evec(SV* obj, double xx, double yy, double zz) {
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->evec[0] = xx;
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->evec[1] = yy;
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->evec[2] = zz;
}

void _xivec(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < 3; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->xivec[i])));
  }
  Inline_Stack_Done;
}
void _set_xivec(SV* obj, double xx, double yy, double zz) {
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->xivec[0] = xx;
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->xivec[1] = yy;
  (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->xivec[2] = zz;
}


/* ----- various things returned by a feff path calculation -- none have setters */
double _edge(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->edge;
}
double _gam_ch(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->gam_ch;
}
double _kf(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->kf;
}
double _mu(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->mu;
}
double _rnorman(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->rnorman;
}
double _rs_int(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->rs_int;
}
double _vint(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->vint;
}
char* _exch(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->exch;
}
char* _version(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->version;
}


double _reff(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->reff;
}
void _ri_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ri[i])));
  }
  Inline_Stack_Done;
}
void _beta_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->beta[i])));
  }
  Inline_Stack_Done;
}
void _eta_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->nleg; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->eta[i])));
  }
  Inline_Stack_Done;
}
void _iz_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < nphx; i++) {
    Inline_Stack_Push(sv_2mortal(newSViv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->iz[i])));
  }
  Inline_Stack_Done;
}


/* ----- number of energy points, may need a setter in the future */
int _ne(SV* obj) {
       return (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne;
}

void _k_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->k[i])));
  }
  Inline_Stack_Done;
}

void _real_phc_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->real_phc[i])));
  }
  Inline_Stack_Done;
}

void _mag_feff_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->mag_feff[i])));
  }
  Inline_Stack_Done;
}

void _pha_feff_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->pha_feff[i])));
  }
  Inline_Stack_Done;
}

void _red_fact_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->red_fact[i])));
  }
  Inline_Stack_Done;
}

void _lam_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->lam[i])));
  }
  Inline_Stack_Done;
}

void _rep_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->ne; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPATH*, SvIV(SvRV(obj))))->rep[i])));
  }
  Inline_Stack_Done;
}
