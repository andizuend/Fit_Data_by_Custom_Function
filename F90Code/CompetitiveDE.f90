!*********************************************************************************************
!*                                                                                           *
!*   Module as an add-on to Differential Evolution (DE) algorithm for global optimization.   *
!*   This module introduces competitive parameter setting into DE according to               *
!*   J. Tvrdik, (2006) Proceedings of the International Multiconference on Computer Science  *
!*   and Information Technology pp. 207–213, ISSN 1896-7094 .                                *
!*                                                                                           *
!*                                                                                           *
!*             (c) Andi Zuend, IACETH, ETH Zurich, 2007 - 2009;                              *
!*       Div. Chem. Engineering, California Institute of Technology, 2009 - 2013             *
!*   Dept. Atmospheric & Oceanic Sciences, McGill University, 2013 --                        *
!*                                                                                           *
!*********************************************************************************************
MODULE CompetitiveDE

USE Mod_NumPrec, ONLY : wp

IMPLICIT NONE
!..
INTEGER,PARAMETER,PRIVATE :: dimF = 5, dimC = 5
INTEGER,PARAMETER,PUBLIC :: dimH = dimF*dimC
REAL(wp),DIMENSION(dimH,2),PUBLIC :: FCcombin
REAL(wp),PARAMETER,PRIVATE :: n0 = 2.0_wp, delta = 1.0_wp/(8.0_wp*dimH)
REAL(wp),DIMENSION(dimH),PUBLIC :: qh    !the probability to choose a certain FCcombination in the next Generation
!..

CONTAINS

!=================================================================================================

    SUBROUTINE FCcombinations

    IMPLICIT NONE
    !..
    INTEGER :: k, fi
    REAL(wp),DIMENSION(dimF),PARAMETER :: F = [ 0.3_wp, 0.5_wp, 0.8_wp, 0.9_wp, 1.0_wp ]
    REAL(wp),DIMENSION(dimC),PARAMETER :: C = [ 0.1_wp, 0.3_wp, 0.5_wp, 0.8_wp, 1.0_wp ]
    !.................

    !create the array with the different combinations of the DE parameters F and C (F_XC and CR_XC):
    fi = 1
    DO k = 1,dimH,dimC !dimC is the increment, so not the usual +1 increments here
        FCcombin(k:k+(dimC-1),1) = F(fi)
        FCcombin(k:k+(dimC-1),2) = C(1:dimC)
        fi = fi+1
    ENDDO

    END SUBROUTINE FCcombinations
    !=================================================================================================

    
    SUBROUTINE FCprobabilities(nh)
    !this subroutine calculates the probability qh for a specific set of F, C parameters  
    !to be chosen, based on their success rate in DE
    IMPLICIT NONE
    REAL(wp),DIMENSION(dimH),INTENT(INOUT) :: nh  !the number of successful individuals using FCcombination h in the current population according to DE comparison
    !...............
    
    !qh is the probability to chose a certain FCcombination in the next generation
    qh = (nh+n0)/(SUM(nh+n0))

    IF (ANY(qh < delta)) THEN !in this case reset the qh's to the initial probability value
        qh = 1.0_wp/REAL(dimH, KIND=wp)
        nh = 0.0_wp !reset count statistics
    ENDIF

    END SUBROUTINE FCprobabilities
    !=================================================================================================

    
    SUBROUTINE FCpop(NPop,Fval,Cval,combPop)
    !this subroutine sets the F and C parameters based on the qh probability distribution 
    !for all Npop individual points
    IMPLICIT NONE
    INTEGER ::  i, k
    INTEGER,INTENT(IN) :: NPop
    INTEGER,DIMENSION(NPop),INTENT(OUT) :: combPop
    REAL(wp),DIMENSION(NPop),INTENT(OUT) :: Fval, Cval !the individual mutation and crossing factors
    REAL(wp),DIMENSION(NPop) :: randvals
    !......................
    
    !generate a set of uniform random numbers between 0 and 1:
    CALL RANDOM_NUMBER(randvals)
    !combine the random numbers with the probabilities qh to get the weighted parameter population:
    DO k = 1,NPop
        DO i = 1,dimH
            IF (SUM(qh(1:i)) > randvals(k)) THEN !found the i value, thus the F and C combination
                Fval(k) = FCcombin(i,1) 
                Cval(k) = FCcombin(i,2)
                combPop(k) = i
                EXIT
            ENDIF
        ENDDO
    ENDDO

    END SUBROUTINE FCpop
    !=================================================================================================

    
END MODULE CompetitiveDE
!============================================ THE END =============================================