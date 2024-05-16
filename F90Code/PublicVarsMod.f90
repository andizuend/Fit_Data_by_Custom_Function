!** Module containing the public variables for the fit data routines. **
MODULE PublicVarsMod

USE Mod_NumPrec, ONLY : wp

IMPLICIT NONE 
!public variables
INTEGER,PUBLIC :: kdatapoints, fcallcount, nvars
REAL(wp),DIMENSION(:),ALLOCATABLE,PUBLIC :: xin, yin, ycalc, amin, amax

END MODULE PublicVarsMod

!============================================================================