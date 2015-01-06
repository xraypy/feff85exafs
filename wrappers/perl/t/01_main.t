#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 43;
use Cwd;

use Xray::FeffPath;

my $here = cwd;
chdir('t') if (cwd !~ m{t\z});

my $epsilon = 1e-4;

my $path = Xray::FeffPath->new();
ok(ref($path) =~ m{FeffPath},								  "object created ".$path);



my $ret = $path->create_path;
ok($ret == $path,									  "called create_path");


my $str = "../../fortran/phase.pad";
$path->phpad($str);
##print $path->phpad, $/;
##(-e $path->phpad) ? print "ok\n" : print "nope\n";
(my $phpad = $path->phpad) =~ s{\s+\z}{};
ok($path->phpad eq $phpad,                                                                "set and read phpad");
ok($path->wrapper->swig_phpad_get eq $phpad ,                                             "low level phpad");

$path->degen(48);
ok($path->degen == 48,									  "set and read degeneracy");
ok($path->wrapper->swig_degen_get == 48,						  "low level degen");



ok($path->Index == 9999,								  "default index");
$path->Index(4);
ok($path->Index == 4,									  "set and read index");
ok($path->wrapper->swig_index_get == 4,							  "low level index");



ok($path->iorder == 2,									  "default iorder");
ok($path->wrapper->swig_iorder_get == 2,						  "low level iorder");



$path->nnnn(1);
ok($path->nnnn,										  "set and read nnnn (boolean)");



$path->evec([0,0,1]);
ok((($path->evec->[0] == 0) and ($path->evec->[1] == 0) and ($path->evec->[2] == 1)),	  "set and read evec");
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


my $ret = $path->path;
ok($ret == $path,									  "called path");

ok($path->exch eq 'H-L exch',                                                             "exchange string ok");
ok(abs($path->edge + 3.82035) < $epsilon,                                                 "edge ok");
ok(abs($path->gam_ch - 1.72919) < $epsilon,                                               "gam_ch ok");
ok(abs($path->kf - 1.824) < $epsilon,                                                     "kf ok");
ok(abs($path->mu + 3.82035) < $epsilon,                                                   "mu ok");
ok(abs($path->rnorman - 2.63173) < $epsilon,                                              "rnorman ok");
ok(abs($path->rs_int - 1.98947) < $epsilon,                                               "rs_int ok");
ok(abs($path->vint + 16.48140) < $epsilon,                                                "vint ok");




ok( (($path->iz->[0] == 29) and ($path->iz->[1] == 29)),                                  "iz array is correct");


ok(( ( abs($path->ri->[0] - 3.61)    < $epsilon) and
     ( abs($path->ri->[1] - 2.55266) < $epsilon) and
     ( abs($path->ri->[2] - 2.55266) < $epsilon) ),					  "ri set correctly");



ok(( ( abs($path->beta->[0] - 135) < $epsilon) and
     ( abs($path->beta->[1] -  90) < $epsilon) and
     ( abs($path->beta->[2] - 135) < $epsilon) ),					  "beta set correctly");



ok(-e 'f3ff0004.dat',									  "feffNNNN.dat file written");
unlink 'f3ff0004.dat';

$path->clear;
ok($path->Index == 9999,								  "path reset");



$path->atom( 0,     0, -3.610, 9);
ok($path->errorcode == 2,                                                                 "error recognized: ipot too big");
$path->clear;



$path->atom( 0,     0, -3.610, -1);
ok($path->errorcode == 1,                                                                 "error recognized: ipot negative");
$path->clear;




$path->atom( 0,     0, -3.610, 1);
$path->atom( 0,     0, -3.710, 1);
ok($path->errorcode == 4,                                                                 "error recognized: atoms too close");
$path->clear;



$path->atom( 0, 0,  0,     1);
$path->atom( 0, 0, -3.61,  1);
$path->path;
ok($path->errorcode == 1,                                                                 "error recognized: begins with absorber");
$path->clear;



$path->atom( 0, 0, -3.61,  1);
$path->atom( 0, 0,  0,     1);
$path->path;
ok($path->errorcode == 2,                                                                 "error recognized: ends with absorber");
$path->clear;



$path->degen(-48);
$path->atom( 0,     0, -3.61,  1);
$path->atom(-1.805, 0, -1.805, 1);
$path->path;
ok($path->errorcode == 4,                                                                 "error recognized: negative degeneracy");
$path->clear;



$path->Index(40000);
$path->atom( 0,     0, -3.61,  1);
$path->atom(-1.805, 0, -1.805, 1);
$path->path;
ok($path->errorcode == 8,                                                                 "error recognized: bad index");
$path->clear;



$path->elpty(-0.5);
$path->atom( 0,     0, -3.61,  1);
$path->atom(-1.805, 0, -1.805, 1);
$path->path;
ok($path->errorcode == 16,                                                                "error recognized: bad elpty");
$path->clear;



$path->iorder(-1);
$path->atom( 0,     0, -3.61,  1);
$path->atom(-1.805, 0, -1.805, 1);
$path->path;
ok($path->errorcode == 32,                                                                "error recognized: bad iorder");
$path->clear;


$path->phpad("foo.bar");
$path->atom( 0,     0, -3.61,  1);
$path->atom(-1.805, 0, -1.805, 1);
$path->path;
ok($path->errorcode == 64,                                                                "error recognized: bad phpad");
#$path->clear;




$path->clean;
#Xray::FeffPathWrapper::FEFFPATH->DISOWN;
undef $path;
#ok(1, "after destruction");
chdir($here);
## fails here,  SIGSEGV (wait status: 139) or SIGFPE (134)

