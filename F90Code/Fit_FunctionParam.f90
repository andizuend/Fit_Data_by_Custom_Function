!*********************************************************************************
!*                                                                               *
!*   This is a simple example program to fit a custom function (parametrization) * 
!*   to a given set of 1-dimensional x <--> y(x) data pairs.                     *
!*   The fit method used is a self-adaptive Best of Random Differential          *
!*   Evolution (BoRDE) variant by Lin et al. (2010, J. Glob. Optim.,             *
!*   doi:10.1007/s10898-010-9535-7), as modified by Andreas Zuend (introducing   *
!*   self-adaptive feature; see Zuend et al. (2010, Atmos. Chem. Phys.,          *
!*   https://doi.org/10.5194/acp-10-7795-2010). Note that other optimization     *
!*   methods may be more efficient and more suitable for a given type of         *
!*   problem. This program simply demonstrates a use case of the BoRDE method.   *
!*                                                                               *
!*                                                                               *
!*   (c) Andi Zuend (andi.zuend@gmail.com),                                      *
!*   Div. Chem. Engineering, California Institute of Technology, 2010            *
!*   Dept. Atmospheric & Oceanic Sciences, McGill University, 2013 --            *
!*                                                                               *
!*********************************************************************************

PROGRAM Fit_FunctionParam

USE Mod_NumPrec, ONLY : wp
USE PublicVarsMod

IMPLICIT NONE

!declare variables and parameters:
INTEGER :: np, i, allocstat, nfitparam, round, istat
REAL(wp) :: fval, a, b
REAL(wp),DIMENSION(:),ALLOCATABLE :: fitpar, amin_stop, amax_stop
CHARACTER(LEN=50) :: filename, dummy
LOGICAL :: exists
!---------------------------------------------------------------------------------------

WRITE(*,*) ""
WRITE(*,*) "////////////////////////////////////////////////////////"
WRITE(*,*) "//  Example program fit function parameters started.  \\"
WRITE(*,*) "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
WRITE(*,*) ""

!read input data file (stored in the main directory of this program):
filename = "InputData.txt"
INQUIRE(FILE = TRIM(filename), EXIST = exists)
IF (exists) THEN
    OPEN (UNIT = 21, FILE = TRIM(filename), STATUS = "OLD")
    READ(21,*) dummy !read first line (the file header)
    READ(21,*) !read empty line 2
    np = 0
    DO  !read the data to a dummy array first to inquire the number of lines
        np = np+1
        READ(21,*,iostat=istat) a, b
        IF (istat /= 0) THEN
            np = np - 1
            EXIT
        ENDIF
    ENDDO
    CLOSE(21)
    !np is now the number of lines
    IF (np > 0) THEN
        ALLOCATE(xin(np), yin(np), ycalc(np), STAT=allocstat)
    ENDIF
    !now read the input data to the corresponding arrays:
    OPEN (UNIT = 21, FILE = TRIM(filename), STATUS = "OLD")
    READ(21,*) dummy !read first line (the file header)
    READ(21,*) !read empty line 2
    DO i = 1,np
        READ(21,*) xin(i), yin(i)
    ENDDO
    CLOSE(21)
ELSE
    WRITE(*,*),""
    WRITE(*,'(A,2X,A)'),"File does not exist or was not found at the designated location: ", filename
    WRITE(*,*),""
    PAUSE
ENDIF
!---------------------------------------------------------------------------------------

!fit the (noisy) data set with the custom function and the objective function defined in the funBoRDE (subject to custom modifications if desired).
!prepare minimization algorithm:
nfitparam = 4   !number of fit parameters; the number of data points, np, needs to be larger than the number of fit parameters!
IF (nfitparam >= np) THEN
    WRITE(*,*),"PROBLEM:"
    WRITE(*,*),"The number of input data points is too small for a reasonable fit with the chosen function (number of fitparameters), nfitparam: ",nfitparam
    WRITE(*,*),"Consider a different custom function (reduced number of function parameters) or provide more input data."
    WRITE(*,*),""
    PAUSE
ENDIF
kdatapoints = np
fcallcount = 0
nvars = nfitparam
ALLOCATE(fitpar(nvars), amin(nvars), amax(nvars), amin_stop(nvars), amax_stop(nvars), STAT=allocstat)
!initialize fit parameters:
fitpar(1:nvars) = 5.0_wp
fitpar(3) = 2.2_wp
fitpar(4) = 0.07_wp

!initial function parameter bounds:
amin(1:nvars) = 1.0E-3_wp  !-100.0_wp
amax(1:nvars) = 1.0E1_wp
amax(3) = 4.0_wp
amax(4) = 1.0_wp

!overall set parameter bounds, never to be exceeded:
amin_stop(1:nvars) = 1.0E-5_wp      !e.g. to keep parameters positive
amax_stop(1:nvars) = 20.0_wp
amax_stop(3) = 20.0_wp
amax_stop(4) = 1.0_wp

DO round = 1,3
    IF (round > 1) THEN
        !revise function parameter bounds based on first estimate:
        amin(1:nvars) = MIN(-ABS(fitpar(1:nvars))*2.0_wp**(round-1), amin(1:nvars))     
        amax(1:nvars) = MAX(ABS(fitpar(1:nvars))*2.0_wp**(round-1), amax(1:nvars))
        !ensure hard limits are fulfilled:
        amin(1:nvars) = MAX(amin(1:nvars), amin_stop(1:nvars))
        amax(1:nvars) = MIN(amax(1:nvars), amax_stop(1:nvars))
    ENDIF
    !--
    CALL MinBoRDE(nvars, amin, amax, fitpar, fval, round)
    !--
    WRITE(*,'(A,2X,I0)'),"end of fitting round no.: ", round
    WRITE(*,*) ""
ENDDO !lo

!call again with best parameterset to calculate ycalc
CALL funBoRDE(fitpar, fval)

WRITE(*,*),""
WRITE(*,*),"finished fitting with Best of Random DE optimization method."
WRITE(*,*),"total number of objective function calls: ", fcallcount
WRITE(*,*),""
WRITE(*,*),"overall objective function value: ", fval
WRITE(*,*),""

!---------------------------------------------------------------------------------------

!output the function parameters and calculated data:
OPEN (UNIT = 22, FILE = "OutputData.txt", STATUS = "UNKNOWN")
WRITE(22,*) "Output of the estimated parameters and function values from the fit routine."
WRITE(22,*) "----------------------------------------------------------------------------"
WRITE(22,*) ""
WRITE(22,'(A,2X,ES15.6)') "objective function value: ", fval
WRITE(22,*) ""
SELECT CASE(nvars)
CASE(3)
    WRITE(22,'(A)') "function parameters: c1,  c2,  c3"
    WRITE(22,*) ""
    WRITE(22,'(3(ES15.6,2X))') fitpar(1), fitpar(2), fitpar(3) 
CASE(4)
    WRITE(22,'(A)') "function parameters: c1,  c2,  c3, c4"
    WRITE(22,*) ""
    WRITE(22,'(4(ES15.6,2X))') fitpar(1), fitpar(2), fitpar(3), fitpar(4) 
CASE(5)
    WRITE(22,'(A)') "function parameters: c1,  c2,  c3, c4, c5"
    WRITE(22,*) ""
    WRITE(22,'(5(ES15.6,2X))') fitpar(1), fitpar(2), fitpar(3), fitpar(4), fitpar(5)
CASE DEFAULT
    WRITE(22,'(A)') "function parameters: c1,  c2"
    WRITE(22,*) ""
    WRITE(22,'(2(ES15.6,2X))') fitpar(1), fitpar(2)
END SELECT
WRITE(22,*) ""
WRITE(22,'(A)') "input data and calculated values (ycalc):"
WRITE(22,*) "    xin,                yin,             ycalc"
WRITE(22,*) ""
DO i = 1,np
    WRITE(22,'(3(ES15.6,3X))') xin(i), yin(i), ycalc(i)
ENDDO
CLOSE(22)

DEALLOCATE(fitpar, amin, amax, STAT=allocstat)
!---------------------------------------------------------------------------------------

!!make a quick, simple DISLIN plot to an x-window showing 
!!the residuals of calculated minus experimental log10(viscosity) vs. mole frac. of water;
!!ifil = 1
!BLOCK
!    REAL(wp) :: xa, xe, xor, xstep, ya, ye, yor, ystep 
!    REAL(wp),PARAMETER :: vtiny = EPSILON(1.0_wp)
!    REAL(wp),DIMENSION(2) :: axrange
!    REAL(wp),DIMENSION(:),ALLOCATABLE :: xray, yray
!    !......................
!    ALLOCATE(xray(np), yray(np))
!    xray(1:np) = xin(1:np)
!    yray(1:np) = (ycalc(1:np) - yin(1:np))
!    !--
!    !!CALL QPLSCA(xray, yray, np)  !one-line plotting with a Dislin quick scatterplot
!    !--
!    !the following is a more elaborate set of plotting options to the screen:
!    CALL METAFL ('xwin')    !'xwin' or 'pdf'
!    CALL SCRMOD ("norev")      !set background color to white and foreground to black
!    CALL DISINI
!    CALL PAGFLL (255)           !white background
!    CALL SETCLR (0)
!    CALL HWFONT                 !for 'xwin' as metafile format
!    CALL TEXMOD ('ON')
!    CALL NAME ('x values', 'X') !x-axis label
!    CALL NAME ('y$_{calc}$ - y$_{inp}$', 'Y')   !y-axis label
!    CALL TITLIN ('Fit residuals', 2)
!    axrange(1:2) = [MINVAL(xray)*0.98_wp, MAX(MAXVAL(xray)*1.01_wp, 1.0_wp)]
!    CALL SETSCL(axrange, 2, 'X')
!    axrange(1:2) = [MINVAL(yray)*1.04_wp, MAXVAL(yray)*1.04_wp]
!    CALL SETSCL(axrange, 2, 'Y')
!    CALL GRAF (xa, xe, xor, xstep, ya, ye, yor, ystep)
!    CALL TITLE
!    CALL HSYMBL (30)            !symbol size
!    CALL SETCLR (20)            !blue color for symbols / curves
!    CALL MARKER (15)            !set symbol type
!    CALL INCMRK(-1)             !plot only symbols
!    CALL CURVE (xray, yray, np) !plot a curve's data
!    CALL ENDGRF
!    CALL DISFIN
!    !--
!    !plot experimental and calculated data (linear y-axis):
!    CALL METAFL ('xwin')        !'xwin' or 'pdf'
!    CALL SCRMOD ("norev")       !set background color to white and foreground to black
!    CALL DISINI
!    CALL PAGFLL(255)            !background color white
!    CALL SETCLR (0)
!    CALL HWFONT
!    CALL TEXMOD ('ON')
!    CALL NAME ('x values', 'X')         !x-axis label
!    CALL NAME ('y_in, y_calc', 'Y')     !y-axis label
!    CALL TITLIN ('Fitted model (curve) vs. input data', 2)
!    axrange(1:2) = [MINVAL(xray)*0.98_wp, MAX(MAXVAL(xray)*1.01_wp, 1.0_wp)]
!    CALL SETSCL(axrange, 2, 'X')
!    axrange(1:2) = [ MIN( MINVAL(yin(1:np))*0.98_wp, MINVAL(ycalc(1:np))*0.98_wp ), &
!        & MAX( MAXVAL(yin(1:np))*1.02_wp, MAXVAL(ycalc(1:np))*1.02_wp ) ]
!    CALL SETSCL(axrange, 2, 'Y')
!    CALL GRAF (xa, xe, xor, xstep, ya, ye, yor, ystep)
!    CALL TITLE
!    CALL HSYMBL (30)        !symbol size
!    CALL SETCLR (20)        !blue color for symbols / curves
!    CALL MARKER (15)        !set symbol type
!    CALL INCMRK(-1)         !plot only symbols
!    CALL CURVE (xray, yin(1:np), np)    !plot exp. data
!    CALL SETCLR (250)       !red
!    CALL INCMRK(0)          !plot curve without symbols
!    CALL CURVE (xray, ycalc(1:np), np)  !plot calc. fit curve data
!    CALL ENDGRF
!    CALL DISFIN
!    DEALLOCATE(xray, yray)
!END BLOCK


IF (exists .AND. np > 0) THEN
    DEALLOCATE(xin, yin, ycalc, STAT=allocstat)
ENDIF

WRITE(*,'(A)'),"best set of function parameters written to the output file: 'OutputData.txt'"
WRITE(*,*) ""
WRITE(*,*) "////////////////////////////////////"
WRITE(*,*) "// program successfully executed. \\"
WRITE(*,*) "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
WRITE(*,*) ""
PAUSE

END PROGRAM Fit_FunctionParam
!=============================== THE END ================================================