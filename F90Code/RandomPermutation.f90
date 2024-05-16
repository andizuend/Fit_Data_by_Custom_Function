!** Subroutine to calculate an array of random number permutations **
!** used in Differential Evolution Algorithms                      **
!** Subroutine slightly modified by Andi Zuend, Caltech, 2010      **

subroutine RandomPermutation(num, randpermut)

use Mod_NumPrec, only : wp

implicit none
integer, intent(in) :: num
integer :: number, i, j, k
integer, dimension(num),intent(out) :: randpermut
real(wp),dimension(num) :: rand2
!...........................................................

call random_number(rand2)
do i=1,num
 number=1
 do j=1,num
    if (rand2(i) > rand2(j)) then
       number=number+1
    end if
 end do
 do k=1,i-1
    if (rand2(i) <= rand2(k) .AND. rand2(i) >= rand2(k)) then
       number=number+1
    end if
 end do
 randpermut(i)=number
end do

end subroutine RandomPermutation

!=================================================================================================