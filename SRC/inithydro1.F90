!--------------------------------------------------
        subroutine inithydro1
!---------------------------------------------------
! ------- modules
        use storage
        implicit none

#ifdef NVFORTRAN        
        real*8 rand 
#else
        real*4 rand ! hack for intel compiler
#endif        
        integer                 ::  i,j
        real(mykind)            ::  x,y, radius2, r
        real(mykind)            ::  shift
        real(mykind), parameter ::  pi=3.141592653589793238462643383279

! set everything to zero
        do j = 1, ny
           do i = 1, nx
              u1(i,j) = 0
              v1(i,j) = 0
              u2(i,j) = 0
              v2(i,j) = 0
           enddo
        enddo

! for a squared box
           radius2 = (radius)**2

! Case 1 (single bubble)
        if (icond == 1) then 
!                
           write(6,*) "INFO: case 1, bubble radius ", sqrt(radius2)
           do j = 0,ny+1
              do i =0,nx+1
!              
                 u1(i,j) = u0
                 v1(i,j) = v0
                 u2(i,j) = u0
                 v2(i,j) = v0
!            
                 if (((i-nx/2)**2+(j-ny/2)**2).lt.radius2) then
                    rhod1(i,j)= rhoin1
                 else
                    rhod1(i,j)=rhoin2
                 endif
              enddo
           enddo
        endif

! Case 2 (many bubbles)
        if (icond == 2) then 
           write(6,*) "INFO: case 2, random noise ", rhoin2
!                
! modified to take into account the possibility of localized peturbation
           do j = 1, ny
              do i = 1, nx
                 u1(i,j) = u0
                 v1(i,j) = v0
                 u2(i,j) = u0
                 v2(i,j) = v0
              enddo
           enddo
!
           write(6,*) "WARNING: magic number"
           do j = 1,ny
              do i =1,nx
#ifdef PWR              
!                 call random_number ( harvest = r )
!                 rhod1(i,j)=0.6931472d0*(1.d0+0.01d0*(r-0.5d0)*2.d0)
!                  rhod1(i,j)=0.5d0*(1.d0+0.01d0*(r-0.5d0)*2.d0)
#else                 
!                 rhod1(i,j)=0.6931472d0*(1.d0+0.01d0*(rand(0)-0.5d0)*2.d0)
                  rhod1(i,j)=0.5d0*(1.d0+0.01d0*(rand(0)-0.5d0)*2.d0)
#endif
              enddo
           enddo
        endif

! Case 3 (two bubbles colliding)
        if (icond == 3) then
           shift = 4*radius/5
           write(6,*) "INFO: case 3, colliding bubble, radius=", & 
     &                 sqrt(radius2)
           write(6,*) "INFO: case 3, colliding bubble, shift=", & 
     &                 shift
!
! left bubble
           do j = 0,ny+1
              do i = 0,nx/2
                 if (((i-nx/4)**2+(j-ny/2-shift)**2).lt.radius2) then
                    rhod1(i,j)= rhoin1
                    u1(i,j) = u0
                    v1(i,j) = v0
                    u2(i,j) = u0
                    v2(i,j) = v0
!
                 else
                    rhod1(i,j)= rhoin2
                 endif
              enddo
           enddo
!
! right bubble                
           do j = 0,ny+1
              do i = nx/2+1, nx+1
                 if (((i-3*nx/4)**2+(j-ny/2+shift)**2).lt.radius2) then
                    rhod1(i,j)= rhoin1
!
                    u1(i,j) = -u0
                    v1(i,j) = -v0
                    u2(i,j) = -u0
                    v2(i,j) = -v0
!
                 else
                    rhod1(i,j)= rhoin2
                 endif
              enddo
           enddo
        endif
!
! Case 4  1-D Band
        if (icond .eq. 4) then
                write(6,*) "INFO: case 4, Flat interface"
           do j = 1, ny
              do i = 1, nx
                 rhod1(i,j)= rhoin1
                 if      (i.lt.(nx/2-radius)) then 
                    rhod1(i,j)= rhoin2
                 else if (i.gt.(nx/2+radius)) then
                    rhod1(i,j)= rhoin2
                 endif
              enddo
           enddo
        endif
!
! Case 5 Taylor-Green Vortex
        if (icond .eq. 5) then 
           write(6,*) "INFO: case 5, Taylor-Green vortex"
           do j = 1, ny
              y = (real(j,mykind)-0.5d0)/real(ny,mykind)        ! 0<y<1 
              do i = 1, nx
                 x = (real(i,mykind)-0.5d0)/real(nx,mykind)     ! 0<x<1 
                 u1(i,j) =+0.1d0*sin(real(2,mykind)*pi*x) & 
     &                          *cos(real(2,mykind)*pi*y)  
                 v1(i,j) =-0.1d0*cos(real(2,mykind)*pi*x) &
     &                          *sin(real(2,mykind)*pi*y)
                 u2(i,j) = u1(i,j)
                 v2(i,j) = v2(i,j)
                 rhod1(i,j) = 1.d0
              enddo
           enddo
        endif
!        
        if ((icond .gt. 5).or.(icond.lt.1)) then 
           write(6,*) "ERROR: option (still) not supported"
           stop
        endif

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine inithydro1"
#endif
        
        end subroutine inithydro1

