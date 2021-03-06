!----------------------------------------------------------
        subroutine out2d
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit none
        
        integer:: i,j
        do j = 1, ny
           do i = 1, nx
              write(51,99) i,j,rhod1(i,j),u1(i,j),v1(i,j)
           enddo
           write(51,'(bn)')
        enddo

        write(51,'(bn)')
        write(51,'(bn)')

 99     format(2I6,3(1x,e13.6))

#ifdef DEBUG
        write(6,*) "Completed subroutine out2d"
#endif
        end subroutine out2d

!----------------------------------------------------------
        subroutine energy
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

! Energy

        g_eff=gnn+gnnn

        en_pot=0.d0
        en_cin=0.d0

        do j = 1, ny
           do i = 1, nx
              en_pot=en_pot+g_eff*psi(i,j)*psi(i,j)*cs2	!/2.d0
              en_cin=en_cin+rhod1(i,j)*cs2
           enddo
        enddo

        en_pot_1=0.d0
        en_pot_2=0.d0

!        write(6,*)'i valori dei c_2 per la subroutine energy sono'
!        write(6,*) c1_2,c2_2,c4_2,c5_2,c8_2


         do j = 1, ny
          do i = 1, nx

! computing forces to find total potential energy

           en_pot_1=en_pot_1+0.5*gnn*psi(i,j)* &
     & (w(1)*c1_2*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     & w(5)*c2_2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2=en_pot_2+0.5d0*gnnn*psi(i,j)* &
     &  (w1*c1_2*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     &   w2*c2_2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+&
     &   w4*c4_2*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+&
     &   w5*c5_2*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+&
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+&
     &  w8*c8_2*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))


         enddo
        enddo

         write(86,98) istep,en_pot_1,en_pot_2,en_pot_1+en_pot_2,en_cin,&
     &                en_pot,(en_pot_1+en_pot_2)/en_cin

98      format(1I6,1x,6(1x,e13.6)) 

! potential energy without "C_i^2" (Prof. Succi 11 - 03 - 2007)

        en_pot_1_bis=0.d0
        en_pot_2_bis=0.d0

         do j = 1, ny
          do i = 1, nx


           en_pot_1_bis=en_pot_1_bis+0.5*gnn*psi(i,j)*&
     &  (4./9.)*psi(i,j)+&
     & (w(1)*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+&
     &  w(5)*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_bis=en_pot_2_bis+0.5*gnnn*psi(i,j)*&
     &  (w0*psi(i,j)+&
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+&
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+&
     &   w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+&
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+&
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+&
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))


         enddo
        enddo

         write(85,96) istep,en_pot_1_bis,en_pot_2_bis, &
     &                en_pot_1_bis+en_pot_2_bis,en_pot

96      format(1I6,1x,4(1x,e13.6))


!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
! Diagnostica a meno di G1 e G2
!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]


        en_pot_ter=0.d0

!        do j = 1, ny
!         do i = 1, nx
!            en_pot_ter=en_pot_ter+psi(i,j)*psi(i,j)*cs2   !/2.d0
!         enddo
!        enddo



        en_pot_1_ter=0.d0
        en_pot_2_ter=0.d0

        en_pot_1_bulk=0.d0
        en_pot_2_bulk=0.d0


!         write(*,*)'i w valgono,',w(0),w(1),w(5)


         do j = 1, ny
          do i = 1, nx

           en_pot_1_ter=en_pot_1_ter+0.5d0*psi(i,j)* &
     & ((4.d0/9.d0)*psi(i,j)+ &
     & w(1)*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     & w(5)*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_ter=en_pot_2_ter+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+ &
     &   w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+ &
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+ &
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+ &
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))



           en_pot_1_bulk=en_pot_1_bulk+0.5d0*psi(i,j)* &
     & ((4.d0/9.d0)*psi(i,j)+ &
     &  w(1)*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &  w(5)*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_bulk=en_pot_2_bulk+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w2*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w4*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w5*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)+ &
     &            psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w8*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))


         enddo
        enddo

97       format(1I6,1x,5(1x,e13.6))

         en_2_1=0.d0
         en_2_2=0.d0
         en_2_1_bulk=0.d0
         en_2_2_bulk=0.d0

         do j = 1, ny
          do i = 1, nx

          en_2_1=en_2_1+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))


          en_2_2=en_2_2+0.5d0*psi(i,j)* &
     &   (w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+ &
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+ &
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+ &
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))



          en_2_1_bulk=en_2_1_bulk+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w2*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))

          en_2_2_bulk=en_2_2_bulk+0.5d0*psi(i,j)* &
     &  (w4*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w5*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)+ &
     &            psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w8*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))


         enddo
        enddo

         write(81,91) istep,en_pot_1_ter,en_pot_2_bulk, &
     &                en_2_1,en_2_2, &
     &                en_2_1_bulk,en_2_2_bulk,en_cin

         write(89,*) istep,gnn*en_pot_1_ter+gnnn*(en_2_1+en_2_2)


91         format(1I6,1x,7(1x,e13.6))

        

!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

        nliq=0

        do i=1,nx
          do j=1,ny

            if(rhod1(i,j).gt.(0.7d0)) then
              nliq=nliq+1
            endif

          enddo
        enddo 


!  number of bubbles...

 
        write(84,*)istep,nliq

        return
        end

