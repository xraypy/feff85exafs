#!/usr/bin/perl

use strict;
use warnings;
use Test::More tests => 21;
use Cwd;

use Xray::FeffPath;

my $here = cwd;
chdir('t') if (cwd !~ m{t\z});

my $epsilon = 1e-4;

my $path = Xray::FeffPath->new();
ok(ref($path) =~ m{FeffPath},								  "object created ".$path);



my $ret = $path->create_path;
ok($ret == $path,									  "called create_path");



$path->deg(48);
ok($path->deg == 48,									  "set and read degeneracy");
ok($path->wrapper->swig_deg_get == 48,							  "low level deg");



ok($path->Index == 9999,								  "default index");
$path->Index(4);
ok($path->Index == 4,									  "set and read index");
ok($path->wrapper->swig_index_get == 4,							  "low level index");



ok($path->iorder == 2,									  "default iorder");
ok($path->wrapper->swig_iorder_get == 2,						  "low level iorder");



$path->nnnn(1);
ok($path->nnnn,										  "set and read nnnn (boolean)");



$path->evec([0,0,1]);
ok((($path->evec->[0] == 0) and ($path->evec->[1] == 0) and ($path->evec->[2] == 1)),	  "set and read xivec");
$path->evec([0,0,0]);
ok((($path->evec->[0] == 0) and ($path->evec->[1] == 0) and ($path->evec->[2] == 0)),	  "unset evec");



$path->xivec([1,1,0]);
ok((($path->xivec->[0] == 1) and ($path->xivec->[1] == 1) and ($path->xivec->[2] == 0)),  "set and read xivec");
$path->xivec([0,0,0]);
ok((($path->xivec->[0] == 0) and ($path->xivec->[1] == 0) and ($path->xivec->[2] == 0)),  "unset xivec");



$path->atom( 0,     0, -3.610, 1);
ok($path->nleg == 2,									  "added first scatterer");



$path->atom(-1.805, 0, -1.805, 1);
ok($path->nleg == 3,									  "added second scatterer");
ok($path->wrapper->swig_nleg_get == 3,							  "low level nleg");



$ret = $path->path;
ok($ret == $path,									  "called path");



ok(( ( abs($path->ri->[0] - 3.61)    < $epsilon) and
     ( abs($path->ri->[1] - 2.55266) < $epsilon) and
     ( abs($path->ri->[2] - 2.55266) < $epsilon) ),					  "ri set correctly");



ok(( ( abs($path->beta->[0] - 135) < $epsilon) and
     ( abs($path->beta->[1] -  90) < $epsilon) and
     ( abs($path->beta->[2] - 135) < $epsilon) ),					  "beta set correctly");



ok(-e 'f3ff0004.dat',									  "feffNNNN.dat file written");

undef $path;
chdir(cwd);
