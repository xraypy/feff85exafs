#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 7;
use Cwd;

use Xray::FeffPath;


my $epsilon = 1e-3;

my $path = Xray::FeffPath->new();
ok(ref($path) =~ m{FeffPath},                                                             "object created ".$path);

### test make_scatterer errors

$path->Index(10);
my $ret = $path->atom(0, 0, -3.61, -1);
ok($ret == 1,                                                                             "ipot negative");

$path->clear;
ok($path->Index == 9999,                                                                  "clear path works");
ok($path->wrapper->_index == 9999,                                                        "clear path works, wrapper");

$ret = $path->atom(0, 0, -3.61, 9);
ok($ret == 2,                                                                             "ipot too big");


$path->clear;

$ret = $path->atom(1.805, 0, 1.805, 1);
$ret = $path->atom(1.805, 0, 1.905, 1);
ok($ret == 4,                                                                             "atoms too close");

$path->clear;
ok($path->nleg == 0,                                                                      "nleg set correctly after clear");


### test make_path errors
