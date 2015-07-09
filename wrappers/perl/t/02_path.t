#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 12;
use Cwd;

use Xray::Feff;


my $epsilon = 1e-3;

my $path = Xray::Feff::Path->new();
ok(ref($path) =~ m{Feff::Path},                                                           "object created ".$path);

my $ret = $path->atom(0, 0, -3.61, 1);
ok((not $ret),                                                                            "added first leg");
ok($path->nleg == 2,                                                                      "nleg=2");
$ret = $path->atom(-1.805, 0, -1.805, 1);
ok((not $ret),                                                                            "added second leg");
ok($path->nleg == 3,                                                                      "nleg = 3");

$path->phpad('../fortran/phase_orig.pad');
$path->degen(48);
$path->Index(4);
$ret = $path->path;
ok((not $ret),                                                                            "made path");

ok($path->ne == 59,                                                                       "correct length of energy grid");
ok(abs($path->reff-4.3577) < $epsilon,                                                    "reff");

ok( (($path->rat->[0]->[0] < $epsilon) and
     ($path->rat->[0]->[1] < $epsilon) and
     (abs($path->rat->[0]->[2]+3.61) < $epsilon)),                                        "first atom coordinates");


ok( ((abs($path->rat->[1]->[0]+1.805) < $epsilon) and
     ($path->rat->[1]->[1] < $epsilon) and
     (abs($path->rat->[1]->[2]+1.805) < $epsilon)),                                       "second atom coordinates");

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
