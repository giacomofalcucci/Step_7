!----------------------------------------------------------
        subroutine out2d_dbg
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit none

        integer:: i,j
        real(sp):: tempPsi
!
        character*19 file_name
        file_name = 'outdbg.xxxxxxxx.vtk'
        write(file_name(8:15),4000) istep
        open(54,file=file_name,status='unknown')

        write(54,'(A26)')'# vtk DataFile Version 2.0'
        write(54,'(A5)')'Campo'
        write(54,'(A6)')'BINARY'
        write(54,'(A25)')'DATASET STRUCTURED_POINTS'
        write(54,'(A11,I10,A1,I10,A1,I10)') 'DIMENSIONS ',nx,' ',ny,' ',1
        write(54,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',1,' ' &
                                                     ,1,' ' &
                                                     ,1
        write(54,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1
        write(54,'(A10,I10)')'POINT_DATA ',nx*ny*1
! write scalar
        write(54,'(A25)')'SCALARS psi float'
        write(54,'(A20)')'LOOKUP_TABLE default'
        close(54)
!
! then write output (binary)
        open(54,file=file_name,status='old', position='append', &
     &          form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
        do j = 1,ny
           do i = 1,nx
              tempPsi = psi(i,j)
              write(54) tempPsi
           end do
        end do
        close(54)

4000    format(i8.8)


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine out2d_dbg"
#endif
 
        end subroutine out2d_dbg

