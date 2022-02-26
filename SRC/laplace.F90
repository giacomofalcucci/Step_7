!----------------------------------------------------------
        subroutine laplace
!----------------------------------------------------------
! ------- modules 
        use storage
        implicit none

        integer:: i,j 
        real(mykind):: pid, pnid, pmax, pmin 
        real(mykind):: rhomax, rhomin, surf_tens 

        p=0.d0

        do j = 1,ny
           do i = 1,nx
              pid  =rhod1(i,j)*cs2
              pnid =0.5d0*(gnn+gnnn)*cs2*psi(i,j)*psi(i,j)
!0.5d0*0.7365d0/3.d0*gnn*psi(i,j)*psi(i,j)+
!                                                          &
!                                                          0.2635d0/3.d0*gnnn*psi(i,j)*psi(i,j)
!                                                          !*psi(i,j)*psi(i,j)
              p(i,j)=pid+pnid
              if(p(i,j).gt.pmax)then
                 pmax = p(i,j)
              endif
              if(p(i,j).lt.pmin)then
                 pmin = p(i,j)
              endif
              write(65,*)i,j,p(i,j)
           enddo
           write(65,'(bn)')
         enddo

         do i=1,nx
            write(66,*)i,p(i,ny/2)
         enddo

         rhoaver=0.d0
         rhomax = -100000000.d0
         rhomin = 100000000.d0
         do j = 1, ny
            do i = 1, nx
               rhomax = max(rhomax,rhod1(i,j))
               rhomin = min(rhomin,rhod1(i,j))
               rhoaver = rhoaver + rhod1(i,j)
            enddo
         enddo

         rhoaver=rhoaver/float(nx*ny)
         write(6,*) 'INFO: mean rho=',rhoaver
         write(6,*) 'INFO: min  rho=',rhomin
         write(6,*) 'INFO: max  rho=',rhomax
         write(6,*) 'INFO: A_1 =',gnn+gnnn
         write(6,*) 'INFO: A_2 =',gnn+1.5d0*gnnn

         radius = ((rhoaver*(nx*ny)-rhomin*(nx*ny))/(3.1415926536d0* &
     &             (rhomax-rhomin)))**(0.5d0)
         surf_tens = (p(nx/2,ny/2)-p(nx/16,ny/16))*radius

         write(6,*) 'INFO: Bubble radius (Case 1)', radius
         write(6,*) 'INFO: Laplace test: gamma = ', surf_tens

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine laplace"
#endif

          end subroutine laplace
 
