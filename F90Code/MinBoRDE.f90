!** the subroutine initialising and using the minimization method. **

SUBROUTINE MinBoRDE(nvars, xl, xu, fitpar, fval, round)

USE Mod_NumPrec, ONLY : wp

IMPLICIT NONE
!...............................
!interface vars:
INTEGER,INTENT(IN) :: nvars, round
REAL(wp),DIMENSION(nvars),INTENT(IN) :: xl, xu
REAL(wp),DIMENSION(nvars),INTENT(INOUT) :: fitpar
REAL(wp),INTENT(OUT) :: fval
!local vars:
INTEGER :: itermax, termcritfit, Npop, nfeval
REAL(wp) :: VTR
!external subroutines to be passed:
EXTERNAL :: funBoRDE
!....................................

!parameters for DE:
SELECT CASE (round)
CASE(1)
    Npop = MIN(MAX(nvars*40, 200), 1000)        !number of individual overall points in each generation
    itermax = 3000      !MIN(100*nvars, 5000)   !maximum number of generations
    termcritfit = 600   !MIN(500, itermax-10)
CASE(2)
    Npop = MIN(MAX(nvars*30, 200), 800) 
    itermax = 1000      !MIN(100*nvars, 5000)
    termcritfit = 400   !MIN(500, itermax-10)
CASE DEFAULT
    Npop = MIN(MAX(nvars*20, 100), 400)         
    itermax = 800       !MIN(100*nvars, 5000)
    termcritfit = 300   !MIN(500, itermax-10)
END SELECT
VTR = -1.0E12_wp

!call the self-adapting "best-of-random" Differential Evolution (BoRDE) algorithm:
!--
CALL BestOfRandomDE(funBoRDE, nvars, xl, xu, VTR, NPop, itermax, fitpar, fval, nfeval, termcritfit) 
!--
!!use the best set of fitpar found:
!CALL funBoRDE(fitpar, fval)     !fval is the objective function output for the best (local) minimum

END SUBROUTINE MinBoRDE
!==========================================================================================