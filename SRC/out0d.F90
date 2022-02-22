!----------------------------------------------------------
        subroutine out0d
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit none

        integer     :: i,j
        integer     :: liquid,vapor, dumb
        real(mykind):: tot_e, tot_u, tot_v
        real(mykind):: densit, paramord
        real(mykind):: rho0, rhomax, rhomin
        real(mykind):: cutoff_l, cutoff_v

        densit = 0.0d0
        tot_u = 0
        tot_v = 0
        tot_e = 0
        rhomin = +1000.d0
        rhomax = -1000.d0
        cutoff_l = rhoin1*0.65           ! liquid cutoff
        cutoff_v = rhoin2*1.35           ! vapor cutoff
        liquid = 0
        vapor = 0

#ifdef NOSHIFT
        fix = zero
#else
        fix = uno
#endif

        do j= 1, ny
           do i = 1, nx
              rho0 = (f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j))     &
     &              +(f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j))     &
     &               +f0(i,j) + fix
!
! computing max/min
              rhomax = max(rho0,rhomax)
              rhomin = min(rho0,rhomin)
!
! computing mean values
              densit = densit + rho0
              tot_u  = tot_u + u1(i,j)
              tot_v  = tot_v + v1(i,j)
              tot_e  = tot_e + ((u1(i,j)*u1(i,j))+(v1(i,j)*v1(i,j)))
!
! computing liqud/vapor
              if (rho0.GT.cutoff_l) then 
                 liquid = liquid + 1 
              else if (rho0.LT.cutoff_v) then
                 vapor = vapor + 1 
              endif
           enddo
        enddo
        densit = densit / dfloat(nx*ny) 
        tot_u = tot_u / dfloat(nx*ny) 
        tot_v = tot_v / dfloat(nx*ny) 
        tot_e = tot_e / dfloat(nx*ny) 
        dumb = nx*ny-vapor-liquid

        write(69,'(I6,4(1x,e13.6))') istep, tot_u, tot_v, tot_e, densit
        write(68,'(I6,2(1x,e13.6),I9,I9,I9)')  &
     &                       istep, rhomax, rhomin, vapor,liquid,dumb



#ifdef DEBUG
        write(6,*) "Completed subroutine out0d"
#endif

        end subroutine out0d

