      program carbonatecalcium
        implicit none
        real(kind=8) :: co2tot=0.0, so4=0.0, no3=0.0, ca, na, pco2
        real(kind=8) :: h=0.0, hco3=0.0, co3=0.0, tol=0.0
        integer :: nca
        integer :: i
        real(kind=8), parameter :: khenry = 29.76
        pco2 = 1.0e-8
        co2tot = pco2 / khenry

        
        tol = 1e-25
        nca = 4
        ca = 0.00
        na = 0.0
        so4 = 0.0
        no3 = 0.0
        
        call carbonate_ca(co2tot, so4, no3,h, hco3,co3,ca, na, tol,nca)
        write(*,*) "Initial CO2 ", co2tot
        write(*,*) "Concentration SO4 ", so4
        write(*,*) "Concentration NO3 ", no3
        write(*,*) "Concentration H+ ", h
        write(*,*) "Concentration CO3 ", co3
        write(*,*) "Concentration HCO3 ", hco3
        write(*,*) "Concentration Ca ", ca
        write(*,*) "pH ", -log10(h)

        open(unit=1, file='concentrations_ca.dat')
        do i=1, 2000
          co2tot = co2tot + 1.0e-5
        call carbonate_ca(co2tot, so4, no3,h, hco3,co3,ca, na, tol,nca)
          write(1,11) co2tot, hco3, co3, ca, na, h, -log10(h)
        end do
11    format(30e12.5,30e12.5,30e12.5,30e12.5,30e12.5,30e12.5,30e12.5)
      end program carbonatecalcium

      subroutine carbonate_ca(co2,so4,no3,h,hco3,co3,ca,na,tol,nca)
      implicit none

      integer :: i , j 
      integer :: nca
      real(kind=8) :: fca_4, fpca_4
      external fca_4, fpca_4

      real(kind=8), parameter :: smallnum = 1.0e-30
      real(kind=8), parameter :: kw = 1.0e-14
      real(kind=8), parameter :: kh = 1.70e-3    !K H2CO3
      real(kind=8), parameter :: k1 = 2.5e-4     !K HCO3
      real(kind=8), parameter :: k2 = 4.69e-11   !K CO3
      real(kind=8), parameter :: ka = 3.2e-09    !K CaCO3
      
      real(kind=8) :: coeffs(nca)
      real(kind=8) :: roots(nca)   !hardcoded changed from n
      real(kind=8) :: so4
      real(kind=8) :: no3
      real(kind=8) :: ca
      real(kind=8) :: na
      real(kind=8) :: tol
      real(kind=8) :: x
      real(kind=8) :: h 
      real(kind=8) :: h2co3
      real(kind=8) :: hco3
      real(kind=8) :: co3
      real(kind=8) :: co2
      logical :: foundroot

      coeffs(1) = - ((kh*k1*k2*co2)**2.)/(ka)
      coeffs(2) = - (kh*k1*k2*co2/(2.*ka))*(kw + k1*kh*co2)
      coeffs(3) = - ((k1*k2*kh*co2)/(2.*ka))*(no3 + na - 2.*so4)
      coeffs(4) = + k2*k1*kh*co2/(2.*ka)

      x = 3.16e-9
      j = 1
      foundroot = .false.

      !call newton(fca_4, fpca_4, x, coeffs, roots, tol, nca)
      call companion_matrix(coeffs, roots, nca)
      
      do i=1,nca
        if (roots(i) > 0.0) then
            foundroot = .true.
            h = roots(i)
        end if 
      end do 

      if(foundroot .eqv. .false.) then
            write(*,*) "Failed to find a positive root"
            write(*,*) roots
            stop
        end if

      if (-log10(h) < 1. .or. -log10(h) > 14.) then
            stop "H+ range out of domain for pH between 1 and 14"
      end if

      if(co2 < smallnum) then
            co2 = 1.0e-10
            write(*,*) "Warning: CO2TOT zero, fixing value at 1.0e-10"
      end if

      
      h2co3 = kh * co2
      hco3 = k1 * h2co3 / h
      co3 = k2 * hco3 / h
      ca = ka / co3 ! if co2tot is very low - lot of ca in solution
      
      end
     
     
      real(kind=8) function fca_4(x, coeffs)  
      implicit none
      real(kind=8) :: x
      real(kind=8) :: coeffs(3)
      !form of the mass balance equation
      fca_4 = x**4 + coeffs(4) * x**3 + coeffs(3) * x**2 +
     +       coeffs(2) * x + coeffs(1)
      return
      end 

      real(kind=8) function fpca_4(x, coeffs) 
      implicit none
      real(kind=8) :: x
      real(kind=8) :: coeffs(3)
      !form of the derivative of the mass balance equation
      fpca_4=4*x**3 + 3*coeffs(4)*x**2 + 2*coeffs(3)* x + coeffs(2)
      return
      end 
      subroutine newton(fca, fpca, x, coeffs, roots, tol, nca)
      implicit none
      real(kind=8) :: fca, fpca
      integer :: nca, j, i
      real(kind=8) :: x
      real(kind=8) :: coeffs(nca)
      real(kind=8) :: roots(nca)
      real(kind=8) :: denominator
      real(kind=8) :: newtonx
      real(kind=8) :: tol
      real(kind=8), parameter :: delta = 0.1
      real(kind=8), parameter :: smallnum = 1.0e-30
      integer, parameter :: MAX_ITER = 300
      logical :: flag
        
      
            do j=1, nca

            do i = 1, MAX_ITER

                  denominator = fpca(x, coeffs)

                  if(abs(denominator) < smallnum) then
                        stop "Error: denominator too small"
                  end if

                  newtonx = x - fca(x, coeffs)/denominator

                  if(abs(newtonx - x) < tol) then
                        flag = .true.
                        x = newtonx + delta
                        roots(j) = newtonx
                        exit
                  else
                  x = newtonx
                  end if
            end do
            
        if(flag .eqv. .false.)  then
                  write(*,*) "Failure to find root. Last: ", newtonx
            end if
      end do

      end subroutine newton
            subroutine companion_matrix(coeffs, roots, nca)
      implicit none

      integer :: i=0
      integer:: nca
      integer :: lwork, info
      real(kind=8) :: coeffs(nca)
      real(kind=8) :: roots(nca)
      real(kind=8) :: matrix(nca, nca) 
      real(kind=8) :: VL(nca, nca)
      real(kind=8) :: VR(nca, nca)
      real(kind=8) :: thing(nca)
      real(kind=8) :: work(nca*3)

      matrix = 0
    
        !companion matrix
      do i=1,nca
        matrix(i,nca) =  -coeffs(i)
        matrix(i+1, i) = 1.0
      end do

      lwork = -1 !hardcoded
      CALL DGEEV( 'N', 'N', nca, matrix, nca, roots, thing, 
     +               VL, nca, VR, nca, work, -1, info )
        lwork = int(work(1))
                
      CALL DGEEV( 'N', 'N', nca, matrix, nca, roots, thing, 
     +                VL, nca, VR, nca, work, lwork, info )
        

      end subroutine companion_matrix