#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 32;
use Cwd;

use Xray::Feff;


my $epsilon = 1e-4;

my $phases = Xray::Feff::Phases->new();
ok(ref($phases) =~ m{Feff::Phases},                                                           "object created ".$phases);

$phases -> jsonfile("../fortran/libpotph.json");
$phases -> _json_pp;

$phases -> phpad("foo.pad");
ok($phases->wrapper->_phpad eq "foo.pad",                                                     "set phpad");

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

ok((($phases->potlbl->[0] =~ m{\ACu\s*\z}) and ($phases->potlbl->[1] =~ m{\ACu\s*\z})),       "potlbl array");
ok(((($phases->wrapper->_potlbl_array)[0] =~ m{\ACu\s*\z}) and (($phases->wrapper->_potlbl_array)[1] =~ m{\ACu\s*\z})),
                                                                                              "potlbl array, wrapper");

# #print join(",", $phases->wrapper->_iz_array), $/;
# #print join(",", @{$phases->iz}), $/;
ok((($phases->lmaxsc->[0] == 2) and ($phases->lmaxsc->[1] == 2)),                             "lmaxsc array");
ok(((($phases->wrapper->_lmaxsc_array)[0] == 2) and (($phases->wrapper->_lmaxsc_array)[1] == 2)),
                                                                                              "lmaxsc array, wrapper");

ok((($phases->xnatph->[0] == 1) and ($phases->xnatph->[1] == 100)),                           "xnatph array");
ok(((($phases->wrapper->_xnatph_array)[0] == 1) and (($phases->wrapper->_xnatph_array)[1] == 100)),
                                                                                              "xnatph array, wrapper");

ok( ( (abs($phases->folp->[0] - 1.15) < $epsilon) and
      (abs($phases->folp->[1] - 1.15) < $epsilon)      ),                                     "folp array");

ok( ( (abs(($phases->wrapper->_folp_array)[0] - 1.15) < $epsilon) and
      (abs(($phases->wrapper->_folp_array)[1] - 1.15) < $epsilon)     ),                      "folp array, wrapper");


## rat is confusing ... it needs to be flattened and unflattened using
## the the interface provided by Inline::C
##
## the next several tests check a few points to verify that this is
## working correctly and specific atom coordinates can be found correctly
my $i = 1;
ok( ( (abs( ($phases->wrapper->_rat_array)[($i-1)*3+0] ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+1] ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+2] ) < $epsilon) ),                       "atom no. $i");

$i = 2;
ok( ( (abs( ($phases->wrapper->_rat_array)[($i-1)*3+0] - 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+1] - 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+2] ) < $epsilon) ),                       "atom no. $i");

$i = 3;
ok( ( (abs( ($phases->wrapper->_rat_array)[($i-1)*3+0] + 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+1] - 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+2] ) < $epsilon) ),                       "atom no. $i");

$i = 13;
ok( ( (abs( ($phases->wrapper->_rat_array)[($i-1)*3+0] ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+1] + 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+2] + 1.805 ) < $epsilon) ),               "atom no. $i");

$i = 23;
ok( ( (abs( ($phases->wrapper->_rat_array)[($i-1)*3+0] + 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+1] - 3.610 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+2] - 1.805 ) < $epsilon) ),               "atom no. $i");

$i = 177;
ok( ( (abs( ($phases->wrapper->_rat_array)[($i-1)*3+0] + 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+1] + 1.805 ) < $epsilon) and
      (abs( ($phases->wrapper->_rat_array)[($i-1)*3+2] + 7.220 ) < $epsilon) ),               "atom no. $i");

ok( ($phases->iphat->[0] == 0 and
     $phases->iphat->[1] == 1 and
     $phases->iphat->[176] == 1),                                                             "iphat");
ok( (($phases->wrapper->_iphat_array)[0] == 0 and
     ($phases->wrapper->_iphat_array)[1] == 1 and
     ($phases->wrapper->_iphat_array)[176] == 1 ),                                            "iphat, wrapper");


# foreach my $i (1,2,3,13,23,177) {
#   print join("|", $i, ($phases->wrapper->_rat_array)[($i-1)*3+0],
# 	     ($phases->wrapper->_rat_array)[($i-1)*3+1],
# 	     ($phases->wrapper->_rat_array)[($i-1)*3+2] ), $/;
# };

# print join("|", 2, rat(1,2), rat(2,2), rat(3,2), iphat(2)), $/;
# print join("|", 3, rat(1,3), rat(2,3), rat(3,3), iphat(3)), $/;
# print join("|", 13, rat(1,13), rat(2,13), rat(3,13), iphat(13)), $/;
# print join("|", 23, rat(1,23), rat(2,23), rat(3,23), iphat(23)), $/;


#$phases->phases;
unlink("foo.pad");
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
