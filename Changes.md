
# Feff85exafs development

## XX November, 2014 (BR)

 1. Compile to shared object libraries
 2. Fortran entry point
 3. C wrapper with error handling
 4. wrote perl wrapper via SWIG
 5. python wrapper via ctypes
 6. yet more improvements to scons build system
 7. documentation

## 3 September, 2014 (BR)

 1. Further improvements to scons build system
 2. Extensive unit testing, see `t/materials/` folder
 3. Modification of source to use json (via [json-fortran](https://github.com/jacobwilliams/json-fortran)) as the intermediate i/o format
 4. No longer writing various intermediate files related to features of Feff not included in feff85exafs
 5. Squelched many compiler warnings, mostly of the "change of type" and "unused parameter" variety

## 17 June, 2014 (BR)

 1.  write SConstruct files for every directory under `src/`
 2.  explicitly identify types in `HEADERS/dim.h` and `HEADERS\const.h`
 3.  fix compilation errors in `EXCH/mpse.f` (mostly fixing undefined
	 types, file had `IMPLICIT NONE`, then failed to define types of
	 many variables)
 
 
 
