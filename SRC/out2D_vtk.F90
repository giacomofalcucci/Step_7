!----------------------------------------------------------
        subroutine out2d_vtk
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit none

        integer:: i,j
        real(sp):: den1
!
        character*19 file_name
        file_name = 'tec_xy.xxxxxxxx.vtk'
        write(file_name(8:15),4000) istep
        open(52,file=file_name,status='unknown')

        write(52,'(A26)')'# vtk DataFile Version 2.0'
        write(52,'(A5)')'Campo'
        write(52,'(A6)')'BINARY'
        write(52,'(A25)')'DATASET STRUCTURED_POINTS'
        write(52,'(A11,I10,A1,I10,A1,I10)') 'DIMENSIONS ',nx,' ',ny,' ',1
        write(52,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',1,' ' &
                                                     ,1,' ' &
                                                     ,1
        write(52,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1
        write(52,'(A10,I10)')'POINT_DATA ',nx*ny*1
! write scalar
        write(52,'(A25)')'SCALARS density float'
        write(52,'(A20)')'LOOKUP_TABLE default'
        close(52)
!
! then write output (binary)
        open(52,file=file_name,status='old', position='append', &
     &          form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
        do j = 1,ny
           do i = 1,nx
              den1 = rhod1(i,j)
              write(52) den1
           end do
        end do
        close(52)

4000    format(i8.8)


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine out2d_vtk"
#endif
 
        end

