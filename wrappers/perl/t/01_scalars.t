#!/usr/bin/perl -w

use strict;
use warnings;
use Test::More tests => 28;
use Cwd;

use Xray::Feff;


my $epsilon = 1e-4;

my $path = Xray::Feff::Path->new();
ok(ref($path) =~ m{Feff::Path},                                                           "object created ".$path);

ok($path->Index == 9999,                                                                  "default index ".$path->Index);
$path->Index(1);
ok($path->Index == 1,                                                                     "set index");
ok($path->wrapper->_index == 1,                                                           "set index, wrapper");

ok($path->degen == 1.0,                                                                   "default degen ".$path->degen);
$path->degen(48);
ok($path->degen == 48,                                                                    "set degen");
ok($path->wrapper->_degen == 48,                                                          "set degen, wrapper");

ok($path->iorder == 2,                                                                    "default iorder ".$path->iorder);
$path->iorder(3);
ok($path->iorder == 3,                                                                    "set iorder");
ok($path->wrapper->_iorder == 3,                                                          "set iorder, wrapper");

ok((not $path->nnnn),                                                                     "default nnnn ".truefalse($path->nnnn));
$path->nnnn(1);
ok($path->nnnn,                                                                           "set nnnn");
ok($path->wrapper->_nnnn,                                                                 "set nnnn, wrapper");

ok((not $path->xdi),                                                                      "default xdi ".truefalse($path->xdi));
$path->xdi(1);
ok($path->xdi,                                                                            "set xdi");
ok($path->wrapper->_xdi,                                                                  "set xdi, wrapper");

ok((not $path->verbose),                                                                  "default verbose ".truefalse($path->verbose));
$path->verbose(1);
ok($path->verbose,                                                                        "set verbose");
ok($path->wrapper->_verbose,                                                              "set verbose, wrapper");

ok((not $path->ipol),                                                                     "default ipol ".truefalse($path->ipol));
$path->ipol(1);
ok($path->ipol,                                                                           "set ipol");
ok($path->wrapper->_ipol,                                                                 "set ipol, wrapper");

ok($path->elpty == 0,                                                                     "default elpty ".$path->elpty);
$path->elpty(0.5);
ok($path->elpty == 0.5,                                                                   "set elpty");
ok($path->wrapper->_elpty == 0.5,                                                         "set elpty, wrapper");

ok($path->phpad eq 'phase.pad',                                                           "default phpad ".$path->phpad);
$path->phpad('../fortran/phase_orig.pad');
ok($path->phpad eq '../fortran/phase_orig.pad' ,                                           "set phpad");
ok($path->wrapper->_phpad eq '../fortran/phase_orig.pad',                                  "set phpad, wrapper");



undef $path;


sub truefalse {
  my ($val) = @_;
  return "true" if $val;
  return "false";
};
