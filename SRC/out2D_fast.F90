!----------------------------------------------------------
        subroutine out2d_fast
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit none

        integer:: i,j
        real(mykind):: cutoff_l
        
                cutoff_l = rhoin1*0.9           ! liquid cutoff

        do j = 1, ny
           do i = 1, nx
                if (rhod1(i,j).GT.cutoff_l) then
                   write(67,98) i,j
                endif
           enddo
        enddo
        write(67,'(bn)')
        write(67,'(bn)')

 98     format(2I6)

#ifdef DEBUG
        write(6,*) "Completed subroutine out2d_fast"
#endif
 
        end

