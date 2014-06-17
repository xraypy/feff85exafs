
# Modifying build system and tracking down build problems

# 17 June, 2014 (BR)

 1.  write SConstruct files for every directory under `src/`
 2.  explicitly identify types in `HEADERS/dim.h` and `HEADERS\const.h`
 3.  fix compilation errors in `EXCH/mpse.f` (mostly fixing undefined
	 types, file had `IMPLICIT NONE`, then failed to define types of
	 many variables)
 
 
 
