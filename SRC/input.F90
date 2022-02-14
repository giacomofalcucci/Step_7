!--------------------------------------------------
        subroutine input
!---------------------------------------------------

! ------- modules
        use storage
        implicit none
!
        integer:: i
!
#ifdef _OPENMP
        use omp_lib
#endif

#ifdef _OPENMP
        INTEGER:: nthreads, threadid
#endif
        open(5,file='muphase.inp')

        print*,'READING: Size X'
        read(5,*)nx

        print*,'READING: Sixe Y'
        read(5,*)ny

        print*,'READING: Number of steps'
        read(5,*)nsteps

        print*,'READING: Number of steps between printing profile'
        read(5,*)nout

        print*,'READING: Number of steps between performing diagnostics'
        read(5,*)ndiag

        print*,'READING: dt'
        read(5,*)dt

        print*,'READING: viscosity'
        read(5,*)visc

        print*,'READING: Coupling gnn'
        read(5,*)gnn

        print*,'READING: rhopsi'
        read(5,*)rhopsi

        print*,' Applied force  (.TRUE. or .FALSE.) ?'
        read(5,*)iforce 

        read(5,*)rhoin1
        print*,'READING: Initial high density', rhoin1
     
        read(5,*)rhoin2
        print*,'READING: Initial low density', rhoin2
     
        print*,'READING: Initial X velocity component'
        read(5,*)u0

        print*,'READING: Initial Y velocity component'
        read(5,*)v0

        print*,'READING: Final velocity for the Poise force'
        read(5,*)uf

        print*,'READING: Linear obstacle ?'
        read(5,*)iobst

        if (iobst) then
            read(5,*)nobst
            print*,'READING: Length of the obstacle (multiple of 2)', nobst
        endif

        print*,'READING: Initial condition (1-4) ?'
        read(5,*)icond

        if(icond.EQ.4) then
            write(6,*) "WARNING: setting gnn=0, single phase"
            gnn = 0
        endif

        print*,'READING: File for output: 5 chars'
        read(5,'(A)')fileout

        read(5,*)dump
        print*,'READING: read populations dump (0 or 1)', dump

        read(5,*)radius
        print*,'READING: read bubble size (0 or 1)', radius

        close(5)

        open(11,file=fileout//'.prof_j.dat')
        open(12,file=fileout//'.prof_i.dat')
        open(51,file=fileout//'.ruv2d')
        open(111,file='dump_pop',status='unknown',form='unformatted')

        open(112,file='dump_rhod1',status='unknown',form='unformatted')

        open(113,file='dump_u1_v1',status='unknown',form='unformatted')

        print*,'*******************************************************'
        print*,'         Lattice BGK model, 2D with 9 velocities'
        print*,'             multiphse code with Shan-Chen EoS'
        print*,'*******************************************************'
        print*,'  developed and released by Prof. GIACOMO FALCUCCI, PhD'
        print*,'  TEST_7  (G.Amati, CINECA)                           '
        print*,'  OpenACC version                                     '
        print*,'               CFD School @ CINECA'
        print*,'                         2021'
#ifdef _OPENMP
!$omp parallel
        nthreads = OMP_GET_NUM_THREADS()
        threadid = OMP_GET_THREAD_NUM()
        if(threadid.eq.0) then
           print*,' Using OpenMP version with threads = ', nthreads
        endif
!$omp end parallel
#endif
!
#ifdef _OPENACC
        print*,' Using OpenACC version                                '
#endif
        print*,'Number of cells :',nx,'*',ny
        print*,'Nsteps :',nsteps
        print*,'Relaxation frequency :',omega
        print*,'Coupling gnn :',gnn
        print*,'Coupling gnnn :',gnnn
        print*,'Applied force :',iforce
        print*,'Initial velocity:',u0, v0
        if (iobst) then
            print*,' Linear Obstacle with length :',nobst
        endif
        write(6,*)'Initial condition', icond
        write(6,*)'Output file :',fileout
        write(6,*) "INFO: mykind=", mykind, "range  =", range(u0)
        write(6,*) "INFO: mykind=", mykind, "huge   =", huge(u0)
        write(6,*) "INFO: mykind=", mykind, "epsilon=", epsilon(u0)
        write(6,*) "INFO: stkind=", stkind, "range  =", range(f0)
        write(6,*) "INFO: stkind=", stkind, "huge   =", huge(f0)
        write(6,*) "INFO: stkind=", stkind, "epsilon=", epsilon(f0)
        print*,'*******************************************************'
#ifdef PWR
! do othing
#else        
        call sleep(1)
#endif        
!
        if (iforce) then
           write(6,*) & 
     &           "INFO: forcing enabled"
           else
           write(6,*) & 
     &           "ERROR: the code is optimized assuming forcing on"
           stop
        endif
! constants

        cs2  = 1.0d0 / 3.0d0
        cs22 = 2.0d0 * cs2
        cssq = 2.0d0 / 9.0d0

! input weights and storage in w(o:npop-1) array

        w1 =4.d0/21.d0/3.d0
        w2 =4.d0/45.d0/3.d0

        w4 =1.d0/60.d0/3.d0
        w5 =2.d0/315.d0/3.d0
        w8 =1.d0/5040.d0/3.d0

        w0=1.d0-(4.d0*w1+4.d0*w2+4.d0*w4+8.d0*w5+4.d0*w8)

        c1_2=1.d0
        c2_2=2.d0

        c4_2=4.d0
        c5_2=5.d0
        c8_2=8.d0 
        
        w(0) = 4.d0/9.d0
        do i = 1, 4
           w(i) = 1.d0/9.d0
           w(i+4) = 1.d0/36.d0
        end do

! reduced density
        den = rhoin/float(npop) 

! scaling
        dx = dt

! calculation of omega
        omega = 1.d0/(3.*visc*(dt*dt)/(dx*dx) + 0.5*dt)

!	visc = (1.0d0 / omega - 0.5d0) * cs2
        print*,'INFO: Viscosity, omega ',visc,omega

! calculation of the constant applied force
        fpois = 8.0d0 * visc * uf / dfloat(ny) / dfloat(ny)
        fpois = rhoin*fpois/6.  ! # of biased populations
        print*,'INFO: Intensity of the applied force ',fpois

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine input"
#endif
       
        end subroutine input
