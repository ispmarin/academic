module carbonate
    implicit none

    integer, parameter, private :: MAX_ITER = 200
    real(kind=8), parameter, private :: epsilon = 1.0e-30
    real(kind=8), parameter :: kw = 1.0e-14
    real(kind=8), parameter :: kh = 1.70e-3
    real(kind=8), parameter :: k1 = 2.5e-4 ! WARNING: H2CO3 or H2CO3 + CO2 dissolved
    real(kind=8), parameter :: k2 = 4.69e-11

contains
	subroutine set_coeff(coeffs, co2tot, so4, no3)
	    implicit none
	    real(kind=8), dimension(:), intent(out) :: coeffs(:)
	    real(kind=8), intent(in) :: co2tot
	    real(kind=8), intent(in) :: so4
	    real(kind=8), intent(in) :: no3

	    coeffs(1) = - 2 * k2 * k1 * kh * co2tot
	    coeffs(2) = - (kw + k1 * kh * co2tot)
	    coeffs(3) = - (2 * so4 + no3)

	end subroutine set_coeff

	function f(x, coeffs)  result (y)
	    implicit none
	    real(kind=8), intent(in) :: x
	    real(kind=8) :: y
	    real(kind=8), dimension(:), intent(in) :: coeffs(:)
	    !form of the mass balance equation
	    y = x**3 + coeffs(3) * x**2 + coeffs(2) * x + coeffs(1)
	end function f

	function fp(x, coeffs) result (y)
	    implicit none
	    real(kind=8), intent(in) :: x
	    real(kind=8) :: y  
	    real(kind=8), dimension(:), intent(in) :: coeffs(:)
	    !form of the derivative of the mass balance equation
	    y=3*x**2 + 2*(coeffs(3))* x + coeffs(2)
	end function fp


    subroutine newton(f, fp, h, hco3, co3, co2tot, tol, flag, coeffs, n)
    	implicit none
    	
    	interface
    		function f(x, coeffs) result (y)
    			real(kind=8), intent(in) :: x
    			real(kind=8) :: y
    			real(kind=8), dimension(:), intent(in) :: coeffs(:)
    		end function f

    		function fp(x, coeffs) result (y)
    			real(kind=8), intent(in) :: x
    			real(kind=8) :: y
    			real(kind=8), dimension(:), intent(in) :: coeffs(:)
    		end function fp
    	end interface

    	integer :: i = 0, j = 1
        integer, intent(in) :: n
    	real(kind=8), intent(in) :: tol
    	real(kind=8) :: denominator
    	real(kind=8) :: newtonx
        real(kind=8) :: x
    	logical, intent(out) :: flag
    	real(kind=8), dimension(:), intent(in) :: coeffs(:)
        real(kind=8) :: roots(n)
        real(kind=8) :: delta = 1.0
        real(kind=8), intent(out) :: h 
        real(kind=8) :: h2co3
        real(kind=8), intent(out) :: hco3
        real(kind=8), intent(out) :: co3
        real(kind=8), intent(in) :: co2tot
        
        do while (j < n)

        	do i = 1, MAX_ITER
        		denominator = fp(x, coeffs)
        		if(abs(denominator) < epsilon) then

        			stop "Error: denominator too small"
        		end if

        		newtonx = x - f(x, coeffs)/denominator

        		if(abs(newtonx - x) < tol) then
        			flag = .true.
        			x = newtonx + delta
                    roots(j) = newtonx
                    j = j+1
                    exit
        		end if

        		x = newtonx
        	end do

        	if(flag .eqv. .false.)  then
        		write(*,*) "Failure to find root. Last approximate root was ", newtonx
        	end if
        end do

        do i=1,n
            if (roots(i) > 0.0) then
                h = roots(i)!-log10(roots(i))
            end if 
        end do 

        h2co3 = kh * co2tot
        hco3 = k1 * h2co3 / h
        co3 = k2 * hco3 / h

    end subroutine newton



end module

