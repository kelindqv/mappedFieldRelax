# mappedFieldRelax
A laid back version of the mappedField boundary condition for OpenFOAM

Usage: Same as the standard OpenFOAM BC mappedField. Use mappedFieldRelax 
on velocity (updated when (iter % 1000) = 0) and mappedFieldRelaxStag on
pressure (updated 500 iterations after velocity).
