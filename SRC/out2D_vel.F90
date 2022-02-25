!----------------------------------------------------------
        subroutine out2d_vel
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit none

        integer:: i,j
        real(sp):: u,v
!
        character*19 file_name
        file_name = 'outvel.xxxxxxxx.vtk'
        write(file_name(8:15),4000) istep
        open(53,file=file_name,status='unknown')

        write(53,'(A26)')'# vtk DataFile Version 2.0'
        write(53,'(A5)')'Campo'
        write(53,'(A6)')'BINARY'
        write(53,'(A25)')'DATASET STRUCTURED_POINTS'
        write(53,'(A11,I10,A1,I10,A1,I10)') 'DIMENSIONS ',nx,' ',ny,' ',1
        write(53,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',1,' ' &
                                                     ,1,' ' &
                                                     ,1
        write(53,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1
        write(53,'(A10,I10)')'POINT_DATA ',nx*ny*1
! write scalar
        write(53,'(A23)')'VECTORS velocity float'
        close(53)
!
! then write output (binary)
        open(53,file=file_name,status='old', position='append', &
     &          form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
        do j = 1,ny
           do i = 1,nx
              u = u1(i,j)
              v = v1(i,j)
              write(53) u, v, 0         ! keep this format, not 0.0
           end do
        end do
        close(53)

4000    format(i8.8)


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine out2d_vel"
#endif
        end subroutine out2d_vel

