#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 15;
use Cwd;

use Xray::Feff;


my $epsilon = 1e-3;

my $path = Xray::Feff::Path->new();
ok(ref($path) =~ m{Feff::Path},                                                     "object created ".$path);

### test make_scatterer errors

$path->Index(10);
my $ret = $path->atom(0, 0, -3.61, -1);
ok($ret == 1,                                                                       "ipot negative");

$path->clear;
ok($path->Index == 9999,                                                            "clear path works");
ok($path->wrapper->_index == 9999,                                                  "clear path works, wrapper");

$ret = $path->atom(0, 0, -3.61, 9);
ok($ret == 2,                                                                       "ipot too big");


$path->clear;

$ret = $path->atom(1.805, 0, 1.805, 1);
$ret = $path->atom(1.805, 0, 1.905, 1);
ok($ret == 4,                                                                       "atoms too close");

$path->clear;
ok($path->nleg == 0,                                                                "nleg set correctly after clear");


### test make_path errors
$path->phbin('../fortran/phase_orig.pad');
$ret = $path->atom(0, 0, 0, 1);
$ret = $path->atom(1.805, 0, 1.805, 1);
$ret = $path->path;
ok($ret < 64,                                                                       "phbin alias works");
ok($ret == 1,                                                                       "first atom is absorber");
$path->clear;

$path->phpad('../fortran/phase_orig.pad');
$ret = $path->atom(1.805, 0, 1.805, 1);
$ret = $path->atom(0, 0, 0, 1);
$ret = $path->path;
ok($ret == 2,                                                                       "last atom is absorber");
$path->clear;

$path->phpad('../fortran/phase_orig.pad');
$ret = $path->atom(0, 0, 3.61, 1);
$ret = $path->atom(1.805, 0, 1.805, 1);
$path->degen(-12.0);
$ret = $path->path;
ok($ret == 4,                                                                       "degeneracy negative");
$path->clear;

$path->phpad('../fortran/phase_orig.pad');
$ret = $path->atom(0, 0, 3.61, 1);
$ret = $path->atom(1.805, 0, 1.805, 1);
$path->Index(40000);
$ret = $path->path;
ok($ret == 8,                                                                       "Index too big");
$path->clear;

$path->phpad('../fortran/phase_orig.pad');
$ret = $path->atom(0, 0, 3.61, 1);
$ret = $path->atom(1.805, 0, 1.805, 1);
$path->elpty(-0.5);
$ret = $path->path;
ok($ret == 16,                                                                      "elpty negative");
$path->clear;

$path->phpad('../fortran/phase_orig.pad');
$ret = $path->atom(0, 0, 3.61, 1);
$ret = $path->atom(1.805, 0, 1.805, 1);
$path->iorder(-1);
$ret = $path->path;
ok($ret == 32,                                                                      "iorder negative");
$path->clear;

$path->phpad('../fortran/phase_orig.padX');
$ret = $path->atom(0, 0, 3.61, 1);
$ret = $path->atom(1.805, 0, 1.805, 1);
$ret = $path->path;
ok($ret == 64,                                                                      "phpad not readable");
$path->clear;


undef $path;
