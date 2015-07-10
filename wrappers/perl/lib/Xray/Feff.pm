package Xray::Feff;
use strict;
use warnings;

our $VERSION = '1.00'; # Inline::MakeMake uses /^\d.\d\d$/ as the
                       # pattern for the version number -- note the
                       # two digits to the right of the dot

sub import {
  strict->import;
  warnings->import;

  foreach my $p (qw(Phases Path)) {
    next if $INC{"Xray/Feff/$p.pm"};
    #print "Xray/Feff/$p.pm\n";
    require "Xray/Feff/$p.pm";
  };
};

1;
__END__

=head1 NAME

Xray::Feff - A wrapper for Feff functionality

=head1 VERSION

This documentation refers to feff85exafs.

=head1 SYNOPSIS

This imports L<Xray::Feff::Path> and L<Xray::Feff::Phases> into your
program.  It also imports L<strict> and L<warnings> into your program.

  use Xray::Feff;
  my $phases = Xray::Feff::Phases->new;
  my $path   = Xray::Feff::Path->new;

=head1 DESCRIPTION

One module to import providing access to all Feff-related
functionality.  Once loaded, L<Xray::Feff::Path> and
L<Xray::Feff::Phases> are available to your program.

=head1 INSTALLATION

After you have built and installed I<feff85exafs>, do the following:

  perl Makefile.PL
  make
  make test
  sudo make install

That's it!  Note, though, that building this wrapper B<requires> that
the fortran and C code for I<feff85exafs> be completely compiled and
that the resulting libraries (and other files) be successfully
installed.

=head1 CONFIGURATION AND ENVIRONMENT

L<Feff85exafs|https://github.com/xraypy/feff85exafs> must be
installed.  Specifically, building the Xray::Feff package must have
access to the F<libfeffphases> and F<libfeffpath> libraries and the
F<feffphases.h> and F<feffpath.h> header files.  These are usually
found in F</usr/local/lib> and F</usr/local/include> on unix and in
??? on Windows.

=head1 DEPENDENCIES

=over 4

=item *

L<Inline>

=item *

L<Inline::C>

=item *

L<Moose>

=item *

L<MooseX::NonMoose>

=item *

L<MooseX::Aliases>

=back

=head1 BUGS AND LIMITATIONS

Please report problems as issues at the github site:
https://github.com/xraypy/feff85exafs

Patches are welcome.

=head1 AUTHOR

Bruce Ravel (L<http://bruceravel.github.io/home>)

http://bruceravel.github.io/demeter/


=head1 LICENCE AND COPYRIGHT

To the extent possible, the authors have waived all rights granted by
copyright law and related laws for the code and documentation that
make up the Perl Interface to the feffphases and feffpath libraries.
While information about Authorship may be retained in some files for
historical reasons, this work is hereby placed in the Public Domain.
This work is published from: United States.

Note that the C libraries themselves are NOT public domain, nor is the
Fortran source code for Feff that they rely upon.

Author: Bruce Ravel (bravel AT bnl DOT gov).
Created: 8 July, 2015

=cut
