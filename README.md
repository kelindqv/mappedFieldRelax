# mappedFieldRelax
A laid back version of the mappedField boundary condition for OpenFOAM

Usage: Same as the standard OpenFOAM BC mappedField. Use mappedFieldRelax 
on velocity (updated when (iter % period) = 0) and mappedFieldRelaxStag on
pressure (updated when (iter % (period + lag) = 0), i.e. after velocity).

Install with wmake and add
libs ( "libmappedFieldRelax.so" );
to your system/controlDict
