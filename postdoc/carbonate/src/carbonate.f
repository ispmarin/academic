      program carbonatetest
        implicit none
        real(kind=8) :: co2tot=0.0, so4=0.0, no3=0.0, ca, na, pco2
        real(kind=8) :: h=0.0, hco3=0.0, co3=0.0, tol=0.0, h2co3
        integer :: nca
        integer :: i
        real(kind=8), parameter :: khenry = 29.76
        pco2 = 2.5 !3.5e-4
        co2tot = pco2 / khenry


        tol = 1e-25
        nca = 3
        ca = 0.4
        na = 0.1
        so4 = 0.41

        call carbonate(co2tot,so4,no3,h,h2co3,hco3,co3,ca,na,tol,nca)
        write(*,*) "initial CO2 ", co2tot
        write(*,*) "Concentration SO4 ", so4
        write(*,*) "Concentration NO3 ", no3
        write(*,*) "Concentration H+ ", h
        write(*,*) "Concentration H2CO3 ", h2co3
        write(*,*) "Concentration CO3 ", co3
        write(*,*) "Concentration HCO3 ", hco3
        write(*,*) "Concentration Ca ", ca
        write(*,*) "pH ", -log10(h)

        open(unit=1, file='concentrations.dat')
        do i=1, 10000
          co2tot = co2tot + 1.0e-5
c          call carbonate(co2tot,so4,no3,h,h2co3,hco3,co3,ca,na,tol,nca)
c          write(1,11) co2tot, h2co3,hco3, co3, ca, na, h, -log10(h)
        end do
11    format(30e12.5,30e12.5,30e12.5,30e12.5,30e12.5,
     +            30e12.5,30e12.5,30e12.5)
      end program carbonatetest


      subroutine carbonate(co2tot,so4,no3,h,h2co3,hco3,
     +           co3,ca,na,tol,nca)
      implicit none

      integer :: i , j
      integer :: nca
      real(kind=8) :: fca, fpca
      external fca, fpca

      real(kind=8), parameter :: smallnum = 1.0e-30
      real(kind=8), parameter :: kw = 1.0e-14
      real(kind=8), parameter :: kh = 1.70e-3
      real(kind=8), parameter :: k1ca = 2.5e-4 !H2CO3 + CO2 dissolved
      real(kind=8), parameter :: k2ca = 4.69e-11
      real(kind=8), parameter :: ka = 3.2e-09   !K dissociation CaCo3
      real(kind=8) :: coeffs(nca)
      real(kind=8) :: roots(nca)   !hardcoded changed from n
      real(kind=8) :: tol
      real(kind=8) :: x
      real(kind=8) :: h
      real(kind=8) :: so4
      real(kind=8) :: no3
      real(kind=8) :: ca
      real(kind=8) :: na
      real(kind=8) :: h2co3
      real(kind=8) :: hco3
      real(kind=8) :: co3
      real(kind=8) :: co2tot
      real(kind=8) :: acso4
      real(kind=8) :: acca
      real(kind=8) :: acna
      real(kind=8) :: achco3
      real(kind=8) :: acco3
      real(kind=8) :: gammaca
      real(kind=8) :: gammana
      real(kind=8) :: gammaso4
      real(kind=8) :: gammaco3
      real(kind=8) :: gammahco3
      real(kind=8) :: icoeff
      logical :: foundroot

      x = 1e-14
      j = 1
      foundroot = .false.


      icoeff = 0.5*(so4*4 + ca*4 + na + hco3 + co3*4)
      gammaca   = 10**(-0.5085*4*(sqrt(icoeff)/(1.0+sqrt(icoeff)) -
     +  0.3*icoeff))
      gammana   = 10**(-0.5085*1*(sqrt(icoeff)/(1.0+sqrt(icoeff)) -
     +   0.3*icoeff))
      gammaso4  = 10**(-0.5085*4*(sqrt(icoeff)/(1.0+sqrt(icoeff)) -
     +   0.3*icoeff))
      gammahco3 = 10**(-0.5085*1*(sqrt(icoeff)/(1.0+sqrt(icoeff)) -
     +   0.3*icoeff))
      gammaco3  = 10**(-0.5085*4*(sqrt(icoeff)/(1.0+sqrt(icoeff)) -
     +   0.3*icoeff))

c      gammaso4 = 1.0
c      gammana = 1.0
c      gammaca = 1.0
c      gammaco3 = 1.0
c      gammahco3 = 1.0

      acso4 = so4 * gammaso4
      acca = ca * gammaca
      acna = na * gammana
      achco3 = hco3 * gammahco3
      acco3 = co3 * gammaco3


      coeffs(1) = - 2. * k2ca * k1ca * kh * co2tot
      coeffs(2) = - (kw + k1ca * kh * co2tot)
      coeffs(3) = + (-2. * acso4 + 2*acca + acna)
      !coeffs(3) = + (-2. * so4 - no3 + 2*ca + na)



      !call newton(fca, fpca, x, coeffs, roots, tol, nca)
      call companion_matrix(coeffs, roots, nca)

      do i=1,nca
        if (roots(i) > 0.0) then
            foundroot = .true.
            h = roots(i)
        end if
      end do

      if(foundroot .eqv. .false.) then
            stop "Failed to find a positive root"
        end if

      if (-log10(h) < 1. .or. -log10(h) > 14.) then
            stop "H+ range out of domain for pH between 1 and 14"
      end if

      if(co2tot < smallnum) then
            co2tot = 1.0e-10
            write(*,*) "Warning: CO2TOT zero, fixing value at 1.0e-10"
      end if


      h2co3 = kh * co2tot
      hco3 = k1ca * h2co3 / h
      co3 = k2ca * hco3 / h



      end subroutine carbonate



      real(kind=8) function fca(x, coeffs)
      implicit none
      real(kind=8) :: x
      real(kind=8) :: coeffs(3)
      !form of the mass balance equation
      fca = x**3 + coeffs(3) * x**2 + coeffs(2) * x + coeffs(1)
      return
      end

      real(kind=8) function fpca(x, coeffs)
      implicit none
      real(kind=8) :: x
      real(kind=8) :: coeffs(3)
      !form of the derivative of the mass balance equation
      fpca=3*x**2 + 2*(coeffs(3))* x + coeffs(2)
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
      integer, parameter :: MAX_icoeffTER = 300
      logical :: flag


            do j=1, nca

            do i = 1, MAX_icoeffTER

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
      lwork =  int(work(1))

      CALL DGEEV( 'N', 'N', nca, matrix, nca, roots, thing,
     +                VL, nca, VR, nca, work, lwork, info )


      end subroutine companion_matrix