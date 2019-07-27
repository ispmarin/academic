! ============================================================================
! Name        : phreeq-fortran.f90
! Author      : Ivan Marin
! Version     :
! Copyright   : Proprietary
! Description :
! ============================================================================


!**********************************************************************
program phreeqcfortran
    use biogeochem
    implicit none
    include "IPhreeqc.f90.inc"
    interface 
        subroutine advection(t2, uo, um, n_nodes, time_step, temp)
            use biogeochem
            implicit none
            real(kind=8), dimension(:,:) :: t2(:,:)
            real(kind=8), dimension(:,:) :: uo(:,:)
            real(kind=8), dimension(:,:) :: um(:,:)
            real(kind=8), dimension(:,:) :: temp(:,:)
            integer :: n_nodes,  time_step
        end subroutine advection
    end interface


    integer :: id,  max_iter=10,  i, j, n_nodes, solution=0, n_solution, k
    integer :: n_spe, n_ele, n_sub, n_b, n_sol
    character(len=40) :: input_file
    character(len=40) :: phreeqc_database
    character(len=20) filename
    character(len=20) concat
    real(kind=8), allocatable :: t2(:,:)
    real(kind=8), allocatable :: mol_weights(:)
    real(kind=8), allocatable :: um(:,:)
    real(kind=8), allocatable :: uo(:,:)
    real(kind=8), allocatable :: temp(:,:)
    real(kind=8):: time_step = 3
    real(kind=8), allocatable :: species_only_array(:)
    real(kind=8) :: t2temp, uotemp
    integer :: nea
    integer :: kinetics
    integer :: max_comp
    integer :: maxnea
    integer :: max_elem
    integer :: num_selected


    character(len=10) sub(10)
    character(len=10) nameEA(4)
    character(len=10) spe(20)
    character(len=10) xele(20)
    character(len=10) xbm(20)
    character(len=10) xsol(20)
    character(len=10) selected_species(20)

    allocate(mol_weights(20))
    mol_weights = 1
    kinetics = 1
    max_comp=10
    maxnea=10
    max_elem=8


    n_nodes = 1
    n_solution=0
    max_iter = 4
    time_step = 8640
    max_iter = 100

    n_spe = 1
    nea = 1
    n_ele = 0
    n_sub = 1
    n_b = 1
    n_sol = 0
    num_selected = 2

    !xele(1) = "C"
    !xele(2) = "O"
    !xele(3) = "H"

    sub(1) = "Benzene"
    !sub(2) = "Xylene"
    !sub(3) = "Toluene"

    xbm(1) = "Biomass"
    !xbm(2) = "Aerogi"
    !xbm(3) = "Putina"

    nameEA(1) = "O2"
    !spe(1) = "CO2"

    !xsol(1) = "Calcite"
    !xsol(2) = "Aragonite"

    selected_species(1) = "CO3-2"
    selected_species(2) = "O2"

    call set_init_values(n_sub, n_spe, n_ele, &
     n_b, n_sol, sub, nameEA, spe, xele, xbm, xsol, nea, mol_weights,kinetics)


    allocate(t2(max_elem,max_comp))
    allocate(um(max_elem,max_comp))
    allocate(uo(max_elem,maxnea))
    allocate(temp(max_elem,max_comp))
    allocate(species_only_array(num_selected))
    
    call getarg(1, input_file)
    call getarg(2, phreeqc_database)


    if (n_nodes >= max_elem) then
        write(*,*) "not possibru"
        stop
    end if

    do j=1,n_nodes
          call init_phreeqc(id, phreeqc_database, input_file, t2(j,:), um(j,:), uo(j,:))
    end do

    i=0
    solution = 1

    do j=1, n_nodes
        write(concat, '(I12)') j
        filename = "concentration" // trim(adjustl(concat)) // ".data"
        open(unit = j,file = filename, status="replace")

        write(concat, '(I12)') j*10
        filename = "conc_phreeqc" // trim(adjustl(concat)) // ".data"
        open(unit = j*10,file = filename, status="replace")
    end do

    do j=1, n_nodes
        t2temp = t2(j,1) 
        uotemp = uo(j,1)    
        write(j,'(A12, 35e12.4, 35e12.4, 35e12.4)') "0", t2temp, uotemp, um(j,1)
        call print_selected_output(j-1,j*10,"file", 0)
    end do


    !t2(1,1) = 0.001
    !uo(1,2) = 0.10

    !****** main loop ******
    do i=1,max_iter 
        do j=1,n_nodes 
            call bionapltobiogeochem(j-1, solution, time_step,  t2(j,:), um(j,:), uo(j,:))
            k = i * time_step
            t2temp = t2(j,1) 
            uotemp = uo(j,5) 
            write(j,'(I12, 35e12.4, 35e12.4, 35e12.4)') k, t2temp, uotemp, um(j, 1)
        end do
            k = i * time_step
            !call advection(t2, uo, um, n_nodes, k, temp)
    end do

    do j=1,n_nodes
        if (DestroyIPhreeqc(j-1) /= IPQ_OK) then
            call OutputErrorString(j)
            stop
       endif
    end do

    do j=1, n_nodes
        close(j)
    end do


    write(*,*) "I'm done!"
end program phreeqcfortran

subroutine advection(t2, uo, um, n_nodes, time_step, temp)
    use biogeochem
    implicit none
    real(kind=8), dimension(:,:) :: t2(:,:)
    real(kind=8), dimension(:,:) :: temp(:,:)
    real(kind=8), dimension(:,:) :: uo(:,:)
    real(kind=8), dimension(:,:) :: um(:,:)
    integer :: n_nodes
    integer :: j
    integer :: time_step
    real :: delta_t = 10, delta_x = 10, vel = 0.5
    
    !using upwind scheme for 1D
    !CFL condition
    if (abs(vel * delta_t/delta_x) > 1) then
        stop "Upwind scheme unstable"
    end if

    temp(1,:) = t2(1,: ) !fixed concentration on first node
    do j=2, n_nodes
        temp(j,:) = t2(j,:) - (delta_t/delta_x)*vel*(t2(j,:) - t2(j-1,:))
    end do
    
    t2 = temp
    
    do j=1, n_nodes
        write(j,'(I12, 15e22.12, 15e22.12, 15e22.12)') time_step, t2(j,1), uo(j,1), um(j, 1)
!        call selected_output_species(j-1, species_only_array, selected_species, num_selected, time_step)
        !call print_selected_output(j-1,j*10,"file", time_step)
    end do



end subroutine advection
