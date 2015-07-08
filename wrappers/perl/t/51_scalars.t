#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 21;
use Cwd;

use Xray::Feff::Phases;


my $epsilon = 1e-4;

my $phases = Xray::Feff::Phases->new();
ok(ref($phases) =~ m{Feff::Phases},                                                           "object created ".$phases);

$phases -> jsonfile("../fortran/libpotph.json");
$phases -> _json_pp;

ok($phases->nat == 177,                                                                       "nat");
ok($phases->wrapper->_nat == 177,                                                             "nat, wrapper");

ok($phases->nph == 1,                                                                         "nph");
ok($phases->wrapper->_nph == 1,                                                               "nph, wrapper");

ok($phases->ntitle == 1,                                                                      "ntitle");
ok($phases->wrapper->_ntitle == 1,                                                            "ntitle, wrapper");

ok($phases->ihole == 1,                                                                       "ihole");
ok($phases->wrapper->_ihole == 1,                                                             "ihole, wrapper");

ok(abs($phases->gamach - 1.72919) < $epsilon,                                                 "gamach");
ok(abs($phases->wrapper->_gamach - 1.72919) < $epsilon,                                       "gamach, wrapper");

ok(abs($phases->rgrd - 0.05) < $epsilon,                                                      "rgrd");
ok(abs($phases->wrapper->_rgrd - 0.05) < $epsilon,                                            "rgrd, wrapper");

ok((($phases->iz->[0] == 29) and ($phases->iz->[1] == 29)),                                   "iz array");
#print '>>>> ', join("|", $phases->wrapper->_iz_array), $/;
ok(((($phases->wrapper->_iz_array)[0] == 29) and (($phases->wrapper->_iz_array)[1] == 29)),   "iz array, wrapper");

# #print join(",", $phases->wrapper->_iz_array), $/;
# #print join(",", @{$phases->iz}), $/;
ok((($phases->lmaxsc->[0] == 2) and ($phases->lmaxsc->[1] == 2)),                             "lmaxsc array");
ok(((($phases->wrapper->_lmaxsc_array)[0] == 2) and (($phases->wrapper->_lmaxsc_array)[1] == 2)), "lmaxsc array, wrapper");

ok((($phases->xnatph->[0] == 1) and ($phases->xnatph->[1] == 100)),                           "xnatph array");
ok(((($phases->wrapper->_xnatph_array)[0] == 1) and (($phases->wrapper->_xnatph_array)[1] == 100)), "xnatph array, wrapper");

ok( ( (abs($phases->folp->[0] - 1.15) < $epsilon) and
      (abs($phases->folp->[1] - 1.15) < $epsilon)      ),                                     "folp array");

ok( ( (abs(($phases->wrapper->_folp_array)[0] - 1.15) < $epsilon) and
      (abs(($phases->wrapper->_folp_array)[1] - 1.15) < $epsilon)     ),                      "folp array, wrapper");



undef $phases;

#  print join(",", $self->wrapper->_iz_array), $/;
#  print join(",", $self->wrapper->_xnatph_array), $/;

# 29,29
# 1,100
# ntitle  1
# nat  177
# nph  1
# ihole  1
# rscf  -1
# lscf  0
# nscmt  0
# ca  0
# nmix  1
# ecv  -40
# icoul  0
# elpty  0
# ispin  0
# angks  0
# gamach  1.72918818490579
# ixc  0
# vr0  0
# vi0  0
# ixc0  0
# iafolp  0
# rgrd  0.0500000007450581
# iunf  0
# inters  0
# totvol  0
# jumprm  0
# nohole  -1
