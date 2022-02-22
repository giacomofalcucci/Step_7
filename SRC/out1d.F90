!----------------------------------------------------------
        subroutine out1d
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit none

        integer:: i,j

! transverse profiles (j)
        do j = 1,ny
           write(11,99) j,rhod1(nx/2,j),u1(nx/2,j),v1(nx/2,j)
        enddo
        write(11,'(bn)')
        write(11,'(bn)')

! transverse profiles (i)
        do i = 1,nx
           write(12,99) i,rhod1(i,ny/2),u1(i,ny/2),v1(i,ny/2)
        enddo
        write(12,'(bn)')
        write(12,'(bn)')

 99     format(I6,3(1x,e13.6))

#ifdef DEBUG
        write(6,*) "Completed subroutine out1d"
#endif
 
        end

