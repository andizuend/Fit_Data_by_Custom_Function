!** subroutine for the objective function **

SUBROUTINE funBoRDE(fitpar, fvalue)

USE Mod_NumPrec, ONLY : wp
USE PublicVarsMod, ONLY : amin, amax, fcallcount, kdatapoints, nvars, ycalc, yin, xin
USE, INTRINSIC :: IEEE_ARITHMETIC

IMPLICIT NONE

!interface variables:
REAL(wp),DIMENSION(nvars),INTENT(INOUT) :: fitpar
REAL(wp),INTENT(OUT) :: fvalue
!local variables:
INTEGER :: i, k
REAL(wp) :: a, d
REAL(wp),DIMENSION(kdatapoints) :: fvec, weighting
!............................................

!-----------------------------------------------------------
!fit noisy data problem (just an example):
k = kdatapoints                 !number of data rows (points in data set to be fitted)

!set weighting of the datapoints:
weighting(1:k) = 1.0_wp  !change the weightings of the data points -- if applicable 


!-----------------------------------------------------------
!4-param modified sigmoidal function; upper and lower limits fixed to known values. 
IF (fitpar(3) * (1.0_wp/fitpar(1))*fitpar(2) > 300.0_wp) THEN           !deal with potential overflow case
    fitpar(1) = fitpar(1)*2.0_wp
    fitpar(3) = 200.0_wp/( (1.0_wp/fitpar(1))*fitpar(2) )
ENDIF
WHERE (fitpar < amin)
    fitpar = amin
ELSE WHERE (fitpar > amax)
    fitpar = amax
END WHERE
weighting(1:k) = 0.05_wp + [(1.0_wp/REAL(i, KIND=wp), i = 1,k)]         !example with different weightings of data points
a = 1.24_wp
d = 1.08_wp != fitpar(4)
DO i = 1,k !k data points
    !!fitpar(1) = 1.0_wp
    !y = d + (a-d)/(1+(x/c)^b)^g 
    ycalc(i) = d + (a - d)/( 1.0_wp + fitpar(4)*xin(i) + (xin(i)/fitpar(1))**fitpar(2) )**fitpar(3)     !4 param model.
    fvec(i) = ABS( LOG(ycalc(i) / yin(i)) )                                                                 !(ycalc(i)-yin(i))**2.0_wp ! deviation squared
END DO
fvalue = SUM(weighting(1:k)*fvec(1:k))
!-----------------------------------------------------------

!!IF ( IEEE_IS_NAN(fvalue) ) THEN
!!    fvalue = fvalue +1.0E3_wp !use a large value indicating a bad fit parameter vector
!!ENDIF

fcallcount = fcallcount + 1 !count the number of objective function calls. 
  
END SUBROUTINE funBoRDE
!============================================================================