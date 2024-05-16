# Example program for Best-of-Random Differential Evolution (BoRDE) optimization method
This Fortran program provides an example use case to demonstrate the application of the Best of Random Differential Evolution (BoRDE) method with adaptive parameter setting.

The main program entry point is in file Fit_FunctionParam.f90, while the main Differential Evolution algorithm is in file/subroutine BestOfRandomDE, which is called from subroutine MinBoRDE.
The code is not particularly polished or documented, but should provide interested users with an example as a starting point. 

# Origin and some details
The fit method used is a self-adaptive Best of Random Differential Evolution (BoRDE) variant by Lin et al. (2010, J. Glob. Optim., https://doi.org/10.1007/s10898-010-9535-7, with some modifications by Andreas Zuend, who introduced the self-adaptive feature; see information in appendix of Zuend et al. (2010, Atmos. Chem. Phys., https://doi.org/10.5194/acp-10-7795-2010). Differential Evolution (DE) as a derivative-free optimization method was first introduced by Storn, R., and Price, K.V., (1996): "Minimizing the real function of the ICEC'96 contest by differential evolution", IEEE conf. on Evolutionary 
Computation, 842-844. The Fortran 90 DE code used here is based on the version written and distributed by Dr. Feng-Sheng Wang, Department of Chemical Engineering, National Chung Cheng University, Chia-Yi 621, Taiwan. It has only been modified slightly and not fully modernized.
