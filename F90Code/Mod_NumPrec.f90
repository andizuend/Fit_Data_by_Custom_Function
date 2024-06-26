!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module defining compiler-independent numerical precision parameters                *  
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend                                                                         *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
module Mod_NumPrec

implicit none

!define a working precision (wp) level to be used with floating point (real) variables, e.g. 1.0 should be stated as 1.0_wp.
!number_of_digits = desired minimum level of precision in terms of number of floating point decimal digits.
integer,parameter,private :: number_of_digits = 12        
integer,parameter,public  :: wp = selected_real_kind(number_of_digits)

end module Mod_NumPrec