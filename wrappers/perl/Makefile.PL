
BEGIN {
  eval {
    require Inline::MakeMaker;
  };
  if ($@) {
    print "You need to install the 'Inline' perl module before building the Xray::Feff wrappers.
See https://metacpan.org/release/Inline
or install from a package for your system:
   debian/ubuntu:     libinline-perl
   redhat/mandriva:   perl-Inline
   suse:              perl-Inline
   gentoo:            dev-perl/Inline
";
    exit;
  }
}

import Inline::MakeMaker;

WriteMakefile(
	      NAME               => 'Xray::Feff',
	      VERSION_FROM       => 'lib/Xray/Feff.pm', # finds $VERSION
	      PREREQ_PM          => {
				     Inline           => 0.5,
				     Inline::C        => 0.62,
				     Moose            => 2.06,
				     MooseX::NonMoose => 0.22,
				     MooseX::Aliases  => 0.10,
				    },
	      ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
	       (#ABSTRACT_FROM    => 'lib/Xray/Feff/Path.pm', # retrieve abstract from module
		ABSTRACT         => 'Perl wrappers libfeffphases and libfeffpath',
		AUTHOR           => 'Bruce Ravel <bravel@bnl.gov>') : ()),
	      LICENSE            => 'unrestricted',
	      CONFIGURE_REQUIRES =>  {
				      'Inline::MakeMaker'     => 0.45,
				      'ExtUtils::MakeMaker'   => 6.52,
				     },
	      #'test' => {'TESTS' => 't/*.t t/gooddata/*.t t/baddata/*.t t/writer/*.t' },
);

# To the extent possible, the authors have waived all rights granted
# by copyright law and related laws for the code and documentation
# that make up the Perl Interface to the FeffPath library.  While
# information about Authorship may be retained in some files for
# historical reasons, this work is hereby placed in the Public Domain.
# This work is published from: United States.
#
# Note that the feffpath and onepath libraries themselves are NOT public
# domain, nor is the Fortran source code for Feff that it relies upon.
#
# Author: Bruce Ravel (bravel AT bnl DOT gov).
# Last update: 13 February, 2015
