#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 5;
use Cwd;

use Xray::Feff;


my $epsilon = 1e-7;

my $path = Xray::Feff::Path->new();
ok(ref($path) =~ m{Feff::Path},                                                           "object created ".$path);


$path->evec([1,1,0]);
ok( ((abs($path->evec->[0] - 1.0) < $epsilon) and
     (abs($path->evec->[1] - 1.0) < $epsilon) and
     ($path->evec->[2] < $epsilon) ),                                                     "polarization vector");

my @vec = $path->wrapper->_evec;
ok( ((abs($vec[0] - 1.0) < $epsilon) and
     (abs($vec[1] - 1.0) < $epsilon) and
     ($vec[2] < $epsilon) ),                                                              "polarization vector, wrapper");


$path->xivec([0,0,1]);
ok( (($path->xivec->[0] < $epsilon) and
     ($path->xivec->[1] < $epsilon) and
     (abs($path->xivec->[2] - 1.0) < $epsilon)),                                           "ellipticity vector");

@vec = $path->wrapper->_xivec;
ok( (($vec[0] < $epsilon) and
     ($vec[1] < $epsilon) and
     (abs($vec[2] - 1.0) < $epsilon)),                                                    "ellipticity vector, wrapper");

undef $path;
