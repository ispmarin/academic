program carbonatetest
    use carbonate
    use carbonate_lapack

    implicit none
    real(kind=8), parameter :: TOL = 1.0e-20
    real(kind=8) :: x
    logical :: flag

    real(kind=8), allocatable :: coeffs(:)
    real(kind=8), allocatable :: roots(:)
    integer :: n = 3
    real(kind=8) :: co2tot, so4, no3, hco3, co3, h

    co2tot = 1.0e-8
    so4 = 0.0
    no3 = 0.0

    allocate(coeffs(n))
    allocate(roots(n))
    call set_coeff(coeffs, co2tot, so4, no3)
    x = 1.0

    call newton(f, fp, h, hco3, co3, co2tot, tol, flag, coeffs, n)
    write(*,*) "Input CO2 concentration: ", co2tot 
    write(*,*) "Calculated concentrations: H+", h 
    write(*,*) "Calculated concentrations: HCO3-", hco3
    write(*,*) "Calculated concentrations: CO3-2", co3
    write(*,*) "Calculated concentrations: pH", -log10(h)
!    call companion_matrix(coeffs, roots, n)

end program carbonatetest

