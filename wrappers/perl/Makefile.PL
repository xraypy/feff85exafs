use lib '.';
use ExtUtils::MakeMaker;
use Config;
use File::Copy;

# my $feffdirs = join(" ",
# 		    '-L/home/bruce/git/feff85exafs/src/GENFMT',
# 		    '-L/home/bruce/git/feff85exafs/src/COMMON',
# 		    '-L/home/bruce/git/feff85exafs/src/MATH',
# 		    '-L/home/bruce/git/feff85exafs/src/PAR',
# 		   );
my $feffdirs = '-L/usr/local/lib';
my $fefflibs = join(" ", qw(-lfeffpath));
#my $fefflibs = join(" ", qw(-lfeffpath -lonepath -lfeffgenfmt -lfeffcom -lfeffmath -lfeffpar -lgfortran -lm));

WriteMakefile('NAME'      => 'Xray::FeffPathWrapper',
              'SKIP'      => ['test_dynamic'],
              'OBJECT'    => 'src/feffpath_wrap.o',
	      'CCFLAGS'   => join(" ", "-o", "src/feffpath_wrap.o", @Config{qw(ccflags optimize cccdlflags)}, "-I$Config{archlib}/CORE"),
              'LDDLFLAGS' => $Config{lddlflags},
              'LIBS'      =>join(" ", $Config{ldflags}, $fefflibs),
              #'LIBS'      =>join(" ", $feffdirs, $Config{ldflags}, $fefflibs),
	     );


#  gcc -c `perl -MConfig -e 'print join(" ", @Config{qw(ccflags optimize cccdlflags)}, "-I$Config{archlib}/CORE")'` feffpath.c feffpath_wrap.c
#  gcc `perl -MConfig -e 'print $Config{lddlflags}'` feffpath.o feffpath_wrap.o -o FeffPath.so
