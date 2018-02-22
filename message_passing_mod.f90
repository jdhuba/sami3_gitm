module message_passing_mod

    use parameter_mod

! Message passing stuff 
            
! this is a buffer that holds deni(nz,nion), ti(nz,nion), and te(nz) 
 
    real :: tl1s(nz,nf,nion+nion+1), tl1r(nz,nf,nion+nion+1) 
    real :: tr1s(nz,nf,nion+nion+1), tr1r(nz,nf,nion+nion+1) 
 
    integer,parameter :: MASTER = 0 
    logical :: flagit, flagit1, flagit10 
    integer :: taskid, source, dest, numworkers 
 
 


end module message_passing_mod
