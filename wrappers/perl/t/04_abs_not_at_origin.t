#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 9;
use Cwd;

use Xray::Feff;


my $epsilon = 1e-3;

my $path = Xray::Feff::Path->new();
ok(ref($path) =~ m{Feff::Path},                                                           "object created ".$path);
$path->absorber([1.1, 1.1, 1.1]);

my $ret = $path->atom(1.1, 1.1, -2.51, 1);
ok((not $ret),                                                                            "added first leg");
ok($path->nleg == 2,                                                                      "nleg=2");
$ret = $path->atom(-0.705, 1.1, -0.705, 1);
ok((not $ret),                                                                            "added second leg");
ok($path->nleg == 3,                                                                      "nleg = 3");

$path->phpad('../fortran/phase_orig.pad');
$path->degen(48);
$path->Index(4);
$ret = $path->path;
ok((not $ret),                                                                            "made path");

ok($path->ne == 59,                                                                       "correct length of energy grid");

ok( ((abs($path->ri->[0]-3.610) < $epsilon) and
     (abs($path->ri->[1]-2.553) < $epsilon) and
     (abs($path->ri->[2]-2.553) < $epsilon)),                                             "ri array");

ok( ((abs($path->beta->[0]-135) < $epsilon) and
     (abs($path->beta->[1]-90 ) < $epsilon) and
     (abs($path->beta->[2]-135) < $epsilon)),                                             "beta array");

#   0.000    0.000   -3.610  1   3.610  135.000    0.000  29
#  -1.805    0.000   -1.805  1   2.553   90.000    0.000  29
#   0.000    0.000    0.000  0   2.553  135.000    0.000  29

undef $path;
