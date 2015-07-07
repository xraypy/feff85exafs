# Using the C wrappers

The wrappers will be installed when feff85exafs is built.  This folder
contains examples of their use in C programs.

 * [Using the libfeffpath wrapper](README_paths.md)
 * [Using the libfeffphases wrapper](README_phases.md)


## Memory testing

If the wrappers have memory leaks or other memory use problems, the
wrappers for other languages will suffer from problems that are very
hard to understand and track down.  It is, therefore, a *very* good
idea to use a tool like [Valgrind](http://valgrind.org/) to test
`makepath`, `errors`, and `makephases` whenever changes are made to
the libraries.

Using [Valgrind](http://valgrind.org/) like so work wonders:

	~> valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all ./makepath
	~> valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all ./errors
	~> valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all ./makephases

If you can get the programs to run without any error reports from
Valgrind, the wrappers for other languages will have a fighting chance
of working correctly.
