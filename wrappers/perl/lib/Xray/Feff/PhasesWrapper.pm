package Xray::Feff::PhasesWrapper;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
# our @EXPORT_OK = (  );
# our @EXPORT = qw(  );

our $VERSION = '1.00';

use Inline C => 'DATA',
           LIBS => '-lfeffphases',
           VERSION => '1.00',
           NAME => 'Xray::Feff::PhasesWrapper';


## for Windows, will need
##   LIBS => '-LC:\\Program Files\\larch\\dlls\\win32 -LC:\\Program Files\\larch\\dlls\\win64 -lfeffphases'
##   INC  => '-IC:\\Program Files\\larch\\include'
##
## also will need Strawberry with gcc 4.9.  Yipes!

1;
__DATA__



=head1 NAME

Xray::FeffPhasesWrapper - Inline::C interface to libfeffphases

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
make up the Perl Interface to the FeffPhases library.  While information
about Authorship may be retained in some files for historical reasons,
this work is hereby placed in the Public Domain.  This work is
published from: United States.

Note that the feffphases and libpotph libraries themselves are NOT public
domain, nor is the Fortran source code for Feff that it relies upon.

Author: Bruce Ravel (bravel AT bnl DOT gov).
Created: 4 June, 2015

=cut

__C__

#ifndef Newx
#  define Newx(v,n,t) New(0,v,n,t)
#endif

#include <string.h>
#include "feffphases.h"

SV* new(char* class) {
  FEFFPHASES * feffphaseswrapper;
  SV *       obj;
  SV *       obj_ref;

  Newx(feffphaseswrapper, 1, FEFFPHASES);

  obj = newSViv((IV)PTR2IV(feffphaseswrapper));
  obj_ref = newRV_noinc(obj);

  sv_bless(obj_ref, gv_stashpv(class, GV_ADD));
  SvREADONLY_on(obj);
  return obj_ref;
}

int    _natx   (SV* obj) { return natx;   }
int    _nphx   (SV* obj) { return nphx;   }
int    _nheadx (SV* obj) { return nheadx; }
double _bohr   (SV* obj) { return bohr;   }
double _ryd    (SV* obj) { return ryd;    }
double _hart   (SV* obj) { return hart;   }


int _create_phases(SV* obj) {
  int ret;
  ret = create_phases((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj)))));
  return ret;
}
void _clear_phases(SV* obj) {
  clear_phases((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj)))));
}
void _dump(SV* obj) {
  dump_phases((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj)))));
}
void _cleanup(SV* obj) {
  cleanup((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj)))));
}
int _make_phases(SV* obj) {
  int ret;
  ret = make_phases((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj)))));
  return ret;
}
int _polarization_tensor(SV* obj) {
  int ret;
  ret = polarization_tensor((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj)))));
  return ret;
}
int _read_libpotph_json(SV* obj) {
  int ret;
  ret = read_libpotph_json((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj)))));
  return ret;
}



/* error handling */

int _errorcode(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->errorcode;
}
char* _errormessage(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->errormessage;
}


/* scalars */

char* _jsonfile(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->jsonfile;
}
void _set_jsonfile(SV* obj, char* c) {
  strcpy((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->jsonfile, c);
}

char* _phpad(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->phpad;
}
void _set_phpad(SV* obj, char* c) {
  strcpy((INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->phpad, c);
}


int _verbose(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->verbose;
}
void _set_verbose(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->verbose = i;
}


int _ntitle(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ntitle;
}
void _set_ntitle(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ntitle = i;
}

int _nat(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nat;
}
void _set_nat(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nat = i;
}

int _nph(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph;
}
void _set_nph(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph = i;
}


int _ihole(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ihole;
}
void _set_ihole(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ihole = i;
}


double _rscf(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->rscf;
}
void _set_rscf(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->rscf = d;
}


int _lscf(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->lscf;
}
void _set_lscf(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->lscf = i;
}


int _nscmt(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nscmt;
}
void _set_nscmt(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nscmt = i;
}


double _ca(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ca;
}
void _set_ca(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ca = d;
}


int _nmix(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nmix;
}
void _set_nmix(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nmix = i;
}


double _ecv(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ecv;
}
void _set_ecv(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ecv = d;
}


int _icoul(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->icoul;
}
void _set_icoul(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->icoul = i;
}


int _ipol(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ipol;
}
void _set_ipol(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ipol = i;
}


double _elpty(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->elpty;
}
void _set_elpty(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->elpty = d;
}


int _ispin(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ispin;
}
void _set_ispin(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ispin = i;
}


double _angks(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->angks;
}
void _set_angks(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->angks = d;
}


double _gamach(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->gamach;
}
void _set_gamach(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->gamach = d;
}


int _ixc(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ixc;
}
void _set_ixc(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ixc = i;
}


double _vr0(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->vr0;
}
void _set_vr0(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->vr0 = d;
}


double _vi0(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->vi0;
}
void _set_vi0(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->vi0 = d;
}


int _ixc0(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ixc0;
}
void _set_ixc0(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ixc0 = i;
}


int _iafolp(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iafolp;
}
void _set_iafolp(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iafolp = i;
}


double _rgrd(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->rgrd;
}
void _set_rgrd(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->rgrd = d;
}


int _iunf(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iunf;
}
void _set_iunf(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iunf = i;
}


int _inters(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->inters;
}
void _set_inters(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->inters = i;
}


double _totvol(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->totvol;
}
void _set_totvol(SV* obj, double d) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->totvol = d;
}


int _jumprm(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->jumprm;
}
void _set_jumprm(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->jumprm = i;
}


int _nohole(SV* obj) {
  return (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nohole;
}
void _set_nohole(SV* obj, int i) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nohole = i;
}


/* ----- 3-vector valued attributes */

void _evec(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < 3; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->evec[i])));
  }
  Inline_Stack_Done;
}
void _set_evec(SV* obj, double xx, double yy, double zz) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->evec[0] = xx;
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->evec[1] = yy;
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->evec[2] = zz;
}

void _xivec(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < 3; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xivec[i])));
  }
  Inline_Stack_Done;
}
void _set_xivec(SV* obj, double xx, double yy, double zz) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xivec[0] = xx;
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xivec[1] = yy;
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xivec[2] = zz;
}


void _spvec(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < 3; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->spvec[i])));
  }
  Inline_Stack_Done;
}
void _set_spvec(SV* obj, double xx, double yy, double zz) {
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->spvec[0] = xx;
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->spvec[1] = yy;
  (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->spvec[2] = zz;
}

/* potentials arrays getters */

void _potlbl_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVpv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->potlbl[i], 6)));
  }
  Inline_Stack_Done;
}
void _set_potlbl_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    strcpy( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->potlbl[i-1], SvPV(Inline_Stack_Item(i), PL_na) );
  }
  Inline_Stack_Void;
}


void _iz_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iz[i])));
  }
  Inline_Stack_Done;
}
void _set_iz_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iz[i-1] = SvIV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


void _lmaxsc_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->lmaxsc[i])));
  }
  Inline_Stack_Done;
}
void _set_lmaxsc_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->lmaxsc[i-1] = SvIV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


void _lmaxph_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->lmaxph[i])));
  }
  Inline_Stack_Done;
}
void _set_lmaxph_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->lmaxph[i-1] = SvIV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


void _xnatph_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xnatph[i])));
  }
  Inline_Stack_Done;
}
void _set_xnatph_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xnatph[i-1] = SvNV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


void _spinph_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->spinph[i])));
  }
  Inline_Stack_Done;
}
void _set_spinph_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->spinph[i-1] = SvNV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


void _folp_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->folp[i])));
  }
  Inline_Stack_Done;
}
void _set_folp_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->folp[i-1] = SvNV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


void _xion_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nph; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xion[i])));
  }
  Inline_Stack_Done;
}
void _set_xion_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nphx+1) {n = nphx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->xion[i-1] = SvNV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


/* atoms arrays */

void _iphat_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nat; i++) {
    Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iphat[i])));
  }
  Inline_Stack_Done;
}
void _set_iphat_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>natx+1) {n = natx+1;}
  for (i=1; i < n; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->iphat[i-1] = SvIV(Inline_Stack_Item(i));
  }
  Inline_Stack_Void;
}


void _rat_array(SV* obj) {
  int i,j;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i < (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->nat; i++) {
    for (j=0; j < 3; j++) {
      Inline_Stack_Push(sv_2mortal(newSVnv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->rat[i][j])));
    }
  }
  Inline_Stack_Done;
}
void _set_rat_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,j,n;
  n=0;
  j=0;
  for (i=1; i < Inline_Stack_Items; i++) {
    (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->rat[n][j] = SvNV(Inline_Stack_Item(i));
    j = j+1;
    if (j == 3) {
      n = n+1;
      if (n > natx) { break; }
      j = 0;
    }
  }
  Inline_Stack_Void;
}



/* headers */
void _titles_array(SV* obj) {
  int i;
  Inline_Stack_Vars;
  Inline_Stack_Reset;
  for (i=0; i <= (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->ntitle; i++) {
    Inline_Stack_Push(sv_2mortal(newSVpv( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->titles[i], 80)));
  }
  Inline_Stack_Done;
}
void _set_titles_array(SV* obj, ...) {
  Inline_Stack_Vars;
  int i,n;
  n = Inline_Stack_Items;
  if (n>nheadx+1) {n = nheadx+1;}
  for (i=1; i < n; i++) {
    strcpy( (INT2PTR(FEFFPHASES*, SvIV(SvRV(obj))))->titles[i-1], SvPV(Inline_Stack_Item(i), PL_na) );
  }
  Inline_Stack_Void;
}


