!** the custom function to fit the input data: **
FUNCTION funcx(x, n, fitpar)

IMPLICIT NONE

INTEGER(4) :: n
REAL(8),INTENT(IN) :: x
REAL(8),DIMENSION(n),INTENT(IN) :: fitpar
REAL(8) :: funcx, a, b, c, d
!............................................

!assign the function parameters (for simpler reading):
a = fitpar(1)
b = fitpar(2)
c = fitpar(3)
d = fitpar(4)

!calculate the function value funcx (= y(x)) at point x:
funcx = a/(x+b) +c*x + d

END FUNCTION funcx

!============================================================================
