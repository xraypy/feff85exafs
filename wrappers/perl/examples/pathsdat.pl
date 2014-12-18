#!/usr/bin/perl

## This reads a paths.dat file and writes f3ffNNNN.dat files using the
## C interface

## test here:
##   cd to ../t
##   perl -I../blib/lib -I../blib/arch/ ../examples/pathsdat.pl


use strict;
use warnings;
use Xray::FeffPath;

my $path = Xray::FeffPath->new();
$path->nnnn(1);
$path->verbose(1);


open(my $PATHS, '<', 'paths.dat');

my $header = 1;
while (<$PATHS>) {

  if ($_ =~ m{-{4,}}) {
    $header = 0;
    next;
  };
  next if $header;

  if ($_ =~ m{\A\s+(\d+)\s+(\d)\s+(\d+\.\d+)\s+}) {
    $path->Index($1);
    $path->deg($3);
    my $line = <$PATHS>;	# column label line
    foreach my $l (1 .. $2-1) {	# snarf atoms
      my $line = <$PATHS>;
      my @list = split(" ", $line);
      $path->atom($list[0], $list[1], $list[2], $list[3]);
    };
    $line = <$PATHS>;	# absorber
    $path->path;
    $path->clear;
  };

}

close $PATHS;
undef $path;
