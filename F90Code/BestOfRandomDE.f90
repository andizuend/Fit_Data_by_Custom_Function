!================================================================================================ 
!     subroutines implementing the best of random differential evolution algorithm to search 
!     for the global minimum of a given objective function.
!================================================================================================

!.......................................................................
!    
! Best of Random Differential Evolution Variant (BoRDE)
!
!.......................................................................
!  This Fortran 90 program translates from the original MATLAB 
!  version of differential evolution (DE). This FORTRAN 90 code 
!  has been tested on Compaq Visual Fortran v6.1. 
!  Any users new to the DE are encouraged to read the article of Storn and Price. 
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of the 
!    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary 
!    Computation, 842-844.
!
!  This Fortran 90 program written by Dr. Feng-Sheng Wang 
!  Department of Chemical Engineering, National Chung Cheng University, 
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!
!  and modified by Andreas Zuend, California Institute of Technology (Caltech), 2010,
!  to include the mutation operation "best-of-random" by Lin, Qing and Feng (2010), J. Glob. Optim.
!.........................................................................


!================================================================================================
!=================================  DE main Subroutine===========================================
!================================================================================================

SUBROUTINE BestOfRandomDE(obj, Dim_XC, XCmin, XCmax, VTR, NPop, itermax, bestmem_XC, bestval, nfeval, termcrit)
!.......................................................................
!    
! Differential Evolution for Optimal Control Problems
!
!.......................................................................
!  This Fortran 90 program translates from the original MATLAB 
!  version of differential evolution (DE). This FORTRAN 90 code 
!  has been tested on Compaq Visual Fortran v6.1. 
!  Any users new to the DE are encouraged to read the article of Storn and Price. 
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of the 
!    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary 
!    Comutation, 842-844.
!
!  This Fortran 90 program written by Dr. Feng-Sheng Wang 
!  Department of Chemical Engineering, National Chung Cheng University, 
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!.........................................................................
!                obj : The user provided file/function for evaluating the objective function.
!                      subroutine obj(xc,fitness)
!                      where "xc" is the real decision parameter vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the real decision parameters.
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
!                VTR : The expected fitness value to reach.
!               NPop : Population size.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for real decision parameters.
!              CR_XC : Crossover factor for real decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced if
!                      "refresh < 1".
! bestmen_XC(Dim_XC) : The best real decision parameters.
!            bestval : The best objective function.
!             nfeval : The number of function CALL.
!         method(1) = 0, Fixed mutation scaling factors (F_XC)
!                   = 1, Random mutation scaling factors F_XC = [0, 1]
!                   = 2, Random mutation scaling factors F_XC = [-1, 1] 
!         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
!                        in the mutation operation 
!                   = other, fixed combined factor provided by the user 
!         method(3) = 1, Saving results in a data file.
!                   = other, displaying results only.
!.......................................................................

USE Mod_NumPrec, ONLY : wp
USE CompetitiveDE

IMPLICIT NONE

!...........................................................
!interface variables:
EXTERNAL :: obj
INTEGER,INTENT(IN) :: NPop, Dim_XC, itermax, termcrit
INTEGER,INTENT(OUT) :: nfeval 
REAL(wp),INTENT(in) :: VTR
REAL(wp),DIMENSION(Dim_XC),INTENT(IN) :: XCmin, XCmax
REAL(wp),DIMENSION(Dim_XC),INTENT(INOUT) :: bestmem_XC
REAL(wp),INTENT(out) :: bestval 
!local variables:   
INTEGER :: i, ibest, iter, iterchanged, itermin, rseed, rbest, rv1, rv2
INTEGER,DIMENSION(NPop) :: rot, a1, a2, a3, rt, combPop
INTEGER,DIMENSION(2) :: ind 
REAL(wp),PARAMETER :: zerolimit = EPSILON(1.0_wp) !accuracy with respect to zero for comparisons; use machine precision
REAL(wp) :: tempval
REAL(wp),DIMENSION(NPop,Dim_XC) :: pop_XC, bm_XC, mui_XC, mpo_XC, rand_XC, ui_XC, vi_XC
REAL(wp),DIMENSION(NPop) :: val, Fval, Cval !the Fval, Cval: individual mutation and crossing factors (used with CompetitiveDE)
REAL(wp),DIMENSION(dimH) :: nh
REAL(wp),DIMENSION(Dim_XC) :: bestmemit_XC
REAL(wp),DIMENSION(Dim_XC) :: rand_C1
!...........................................................

!initialize the combinations for use in the DE-algorithm
CALL FCcombinations 

!initialize random number generator
!rseed = 25779443
!CALL RANDOM_SEED(rseed) !constant seed for debugging
CALL RANDOM_SEED() !use time from computer as seed

itermin = MIN(Dim_XC*14, itermax/3)

!!-----Initialize a population --------------------------------------------!!
!check bounds for the first (old_bestmem_XC) guess values:
WHERE (bestmem_XC < XCmin)
    bestmem_XC = XCmin 
ELSE WHERE (bestmem_XC > XCmax)
    bestmem_XC = XCmax 
END WHERE
pop_XC(1,:) = bestmem_XC !the best parameters from the previous solution or just some input values
DO i = 2,NPop !set the rest of the population values by a pseudo random number generator
    CALL RANDOM_NUMBER(rand_C1)
    pop_XC(i,:) =  XCmin+rand_C1*(XCmax-XCmin)
ENDDO
!!--------------------------------------------------------------------------!!

!!------Evaluate fitness functions and find the best member-----------------!!
nfeval = 0
ibest = 1
CALL obj(pop_XC(1,:), val(1))  
bestval = val(1)
nfeval = nfeval+1
DO i = 2,NPop
    CALL obj(pop_XC(i,:), val(i))
    nfeval = nfeval+1
    IF (val(i) + zerolimit < bestval) THEN
        ibest = i
        bestval = val(i)
    ENDIF
ENDDO  	 
bestmemit_XC = pop_XC(ibest,:)
bestmem_XC = bestmemit_XC
!!--------------------------------------------------------------------------!!
bm_XC = 0.0_wp
rot = [(i,i = 0,NPop-1)]
iter = 1 
iterchanged = 1

!initialize the count statistics for the probabilities qh:
nh = 0.0_wp
 
!!------Perform evolutionary computation------------------------------------!! 
 DO WHILE (iter <= itermax)

!!------Mutation operation--------------------------------------------------!!
    CALL RandomPermutation(2,ind)
    CALL RandomPermutation(NPop, a1)
    rt = MOD(rot+ind(1), NPop)
    a2 = a1(rt+1)
    rt = MOD(rot+ind(2), NPop)
    a3 = a2(rt+1)
	
!-----  use competitive probabilities to set the scaling factor F_XC and the mutation parameter CR_XC (implemented by Andi Z.) -----!
    CALL FCprobabilities(nh) !this routine sets the qh probabilities
    CALL FCpop(NPop,Fval,Cval,combPop) !this routine sets the values for Fval, Cval and registers the set in combPop

!---- perform mutation -----------------------------------------!
    DO i = 1,NPop
        !choose the best-of-random vector out of the three randomly selected individuals:
        IF (val(a1(i)) < val(a2(i))) THEN
            rbest = a1(i) !rbest = random best
            rv1 = a2(i)
            IF (val(a3(i)) < val(rbest)) THEN
                rbest = a3(i)
                rv2 = a1(i)
            ELSE
                rv2 = a3(i)
            ENDIF
        ELSE
            rbest = a2(i)
            rv1 = a1(i)
            IF (val(a3(i)) < val(rbest)) THEN
                rbest = a3(i)
                rv2 = a2(i)
            ELSE
                rv2 = a3(i)
            ENDIF
        ENDIF
        !mutation vector with the best of the three as basis and the other two used for the vector difference (this enhances diversity)
        vi_XC(i,:) = pop_XC(rbest,:) + Fval(i)*( pop_XC(rv1,:)- pop_XC(rv2,:) )
    ENDDO
  
!!--------------------------------------------------------------------------!!
!!------Crossover operation-------------------------------------------------!!
    CALL RANDOM_NUMBER(rand_XC)
    mui_XC = 0.0_wp
    mpo_XC = 0.0_wp
    !version with CompetitiveDE:
    DO i = 1,NPop
        WHERE (rand_XC(i,:) < Cval(i))
           mui_XC(i,:) = 1.0_wp
        ELSEWHERE
           mpo_XC(i,:) = 1.0_wp
        ENDWHERE
    ENDDO
	ui_XC = pop_XC*mpo_XC + vi_XC*mui_XC
	
!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!
    DO i = 1,NPop
    !!------Confine each of feasible individuals within the lower-upper bounds-------!!
    !        ui_XC(i,:) =  MAX( MIN(ui_XC(i,:), XCmax), XCmin )
        IF (ANY(ui_XC(i,:) < XCmin)) THEN
            WHERE (ui_XC(i,:) < XCmin)
                ui_XC(i,:) = XCmin+(XCmax-XCmin)*ABS(ui_XC(i,:) / (ui_XC(i,:)+XCmin))   !this is better than just setting the value to the lower bound
            ENDWHERE
        ENDIF
        IF (ANY(ui_XC(i,:) > XCmax)) THEN
            WHERE (ui_XC(i,:) > XCmax)
                ui_XC(i,:) = XCmax-(XCmax-XCmin)*ABS(ui_XC(i,:) / (ui_XC(i,:)+XCmax))
            ENDWHERE
        ENDIF
        CALL obj(ui_XC(i,:), tempval)
        nfeval = nfeval+1
        nh(combPop(i)) = nh(combPop(i)) + 1.0_wp                                        !count successful events for the next competitive selection of F and C values based on probabilities 
        IF (tempval < val(i)) THEN
            pop_XC(i,:) = ui_XC(i,:)
            val(i) = tempval
            IF (tempval+zerolimit < bestval) THEN 
                bestval = tempval
                bestmem_XC = ui_XC(i,:)
                iterchanged = iter
            ENDIF
        ENDIF
    ENDDO
    bestmemit_XC = bestmem_XC
    iter = iter +1
    !generate some intermediate output
    IF ((MOD(iter,20) == 0 .AND. iter < 101) .OR. MOD(iter,200) == 0) THEN
        WRITE(*,'(A,I0.3,A,ES12.5)') "In BoRDE at iteration, ",iter,", best val. so far = ", bestval
    ENDIF
    !check for maximum number of iteration with no change in the best function value found. (the first criterion for termination...)
    IF ( ( (iter-iterchanged > termcrit*Dim_XC .AND. iter > itermin)  .OR. iter > itermax ) ) THEN
        WRITE(*,*) "................................................."
        WRITE(*,*) "exiting Best of Random Diff. Evol. optimization: "
        WRITE(*,*) "no further changes in the best match."
        WRITE(*,*) "iter, iterchanged:  ", iter, iterchanged
        WRITE(*,*) "number of function evaluations: ", nfeval
        WRITE(*,*) "best function value: ", bestval
        WRITE(*,*) "................................................."
        WRITE(*,*) ""
        EXIT !exiting the loop
    ENDIF
    IF (bestval  <  VTR) THEN
	   WRITE(*,'(A)') 'WARNING from BestOfRandomDE: The best fitness is smaller than VTR' 
       EXIT
    ENDIF
 ENDDO
!!------end the evolutionary computation------------------------------!!

END SUBROUTINE BestOfRandomDE
!================================================================================================