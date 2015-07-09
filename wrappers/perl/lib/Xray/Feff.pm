package Xray::Feff;
use strict;
use warnings;

our $VERSION = '1.00'; # Inline::MakeMake uses /^\d.\d\d$/ as the pattern for the version number -- note the two digits to the right of the dot

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

One module to import to gain access to all Feff-related functionality.
Once loaded, L<Xray::Feff::Path> and L<Xray::Feff::Phases> are
available to your program.

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

Copyright (c) 2015 Bruce Ravel (L<http://bruceravel.github.io/home>). All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlgpl>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
