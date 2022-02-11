!
        module storage
!
        implicit double precision(a-h,o-z)
!
        integer, parameter:: sp=kind(1.0)
!        integer, parameter:: dp=selected_real_kind(2*precision(1.0_sp))
!
        integer :: nx
        integer :: ny
#ifdef SINGLEP
! current precision
        integer,parameter:: mykind=selected_real_kind(1*precision(1.0_sp))
!
! storage precison
        integer,parameter:: stkind=selected_real_kind(1*precision(1.0_sp))    
#elif MIXEDP
! current precision
        integer,parameter:: mykind=selected_real_kind(2*precision(1.0_sp))
!
! storage precison
        integer,parameter:: stkind=selected_real_kind(1*precision(1.0_sp))    
#else 
! current precision
        integer,parameter:: mykind=selected_real_kind(2*precision(1.0_sp))
!
! storage precison
        integer,parameter:: stkind=selected_real_kind(2*precision(1.0_sp))
#endif
        integer, parameter :: npop = 9
!
        character*5 fileout
        integer::  nrhout
        logical iforce,iobst
!/phys/ 
        real(mykind) u0,v0,uf,fom  
        real(mykind) f_guo
!/constants/
        real(mykind) cs2,cs22,cssq,rhoin,omega,fpois,den,visc
        real(mykind) w0,w1,w2,w4,w5,w8,gnn,gnnn,rhoaver,dinvrho
        real(mykind) rhopsi,dt,dx,dump,c1_2,c2_2,c4_2,c5_2,c8_2 
        real(mykind), parameter :: cte04 = (4.d0/ 9.d0)
        real(mykind), parameter :: cte09 = (1.d0/ 9.d0)
        real(mykind), parameter :: cte36 = (1.d0/36.d0)
!/count/ 
        integer istep,nout,ndiag,nsteps,nobst,icond
!/arrays to allocate/
        real(stkind), dimension (:,:), allocatable :: f0
        real(stkind), dimension (:,:), allocatable :: f1
        real(stkind), dimension (:,:), allocatable :: f2
        real(stkind), dimension (:,:), allocatable :: f3
        real(stkind), dimension (:,:), allocatable :: f4
        real(stkind), dimension (:,:), allocatable :: f5
        real(stkind), dimension (:,:), allocatable :: f6
        real(stkind), dimension (:,:), allocatable :: f7
        real(stkind), dimension (:,:), allocatable :: f8
!
        real(stkind), dimension (:,:), allocatable :: fp0
        real(stkind), dimension (:,:), allocatable :: fp1
        real(stkind), dimension (:,:), allocatable :: fp2
        real(stkind), dimension (:,:), allocatable :: fp3
        real(stkind), dimension (:,:), allocatable :: fp4
        real(stkind), dimension (:,:), allocatable :: fp5
        real(stkind), dimension (:,:), allocatable :: fp6
        real(stkind), dimension (:,:), allocatable :: fp7
        real(stkind), dimension (:,:), allocatable :: fp8
!
        real(stkind), dimension (:,:), allocatable :: u1
        real(stkind), dimension (:,:), allocatable :: v1
        real(stkind), dimension (:,:), allocatable :: u2
        real(stkind), dimension (:,:), allocatable :: v2
        real(stkind), dimension (:,:), allocatable :: psi
        real(stkind), dimension (:,:), allocatable :: rhod1
        real(stkind), dimension (:,:), allocatable :: rhod2
        real(stkind), dimension (:,:), allocatable :: p
        real(stkind), dimension (:,:), allocatable :: param
!       
        integer, dimension (:,:), allocatable ::  iflag

!/arrays/
        real(mykind), dimension (0:npop-1) ::     w
        real(mykind), dimension (0:npop-1) ::     u_ci
!                
!        
        end module  storage

