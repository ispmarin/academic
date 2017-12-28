module carbonate_lapack    
  implicit none

  contains

    subroutine companion_matrix(coeffs, roots, n)
    	implicit none

    	integer :: i=0
    	integer, intent(in) :: n
    	real(kind=8), dimension(:), intent(in) :: coeffs(n)
    	real(kind=8), dimension(:), intent(out) :: roots(n)
		real(kind=8), dimension(:,:) :: matrix(n, n) 
		integer :: lwork, info
		real(kind=8), allocatable :: VL(:, :)
		real(kind=8), allocatable :: VR(:, :)
		real(kind=8), allocatable :: thing(:)
		real :: work(n*3)
		
		allocate(VL(n, n))		
		allocate(VR(n, n))
		allocate(thing(n))
        matrix = 0
		
        !companion matrix
		do i=1,n
			matrix(i,n) =  -coeffs(i)
			matrix(i+1, i) = 1.0
		end do

      	lwork = -1 !hardcoded
      	CALL DGEEV( 'N', 'N', n, matrix, n, roots, thing, VL, n, VR, n, work, -1, info )
      	lwork = 100!int(work(2))
      	      	
      	CALL DGEEV( 'N', 'N', n, matrix, n, roots, thing, VL, n, VR, n, work, lwork, info )
      	write(*,*) "roots", roots
        
        do i=1,n
            if (roots(i) > 0.0) then
                write(*,*) -log10(roots(i))
            end if 
        end do 

	end subroutine companion_matrix

end module carbonate_lapack