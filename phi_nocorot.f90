
! PHI_NOCOROT.F

! take out corotation potential in sami3
! need to set ntime

    program MAIN

    use parameter_mod
         
!    parameter ( ntime = 12 )

    real, dimension(:,:,:), allocatable :: phi,phi_nocorot

    real :: baltp(nzp1,nfp1,nlt)
    real :: phi_corot(nny),preal(nny)
    real :: temp(nnx,nny)

    open ( unit=497, file='baltpu.dat',    form='unformatted' )
    open ( unit=498, file='phiu.dat',      form='unformatted' )
    open ( unit=499, file='phi_nocoru.dat',form='unformatted' )

    print *,'input ntime'
    read *,ntime

    allocate (phi(nnx,nny,ntime),phi_nocorot(nnx,nny,ntime))


! read in baltp and define preal

    read(497) baltp

! fine max of baltp

    baltpmax = -1.e2

    k    = 1
    do j = 1,nf+1
        do i = 1,nz+1
            if ( baltp(i,j,k) > baltpmax ) then
                baltpmax = baltp(i,j,k)
                imax     = i
                jmax     = j
            endif
        enddo
    enddo
           
    print *,'max ',imax,jmax,baltpmax/re

    do j = 1,nf-1
        imid     = (nz + 2)/2
        preal(j) = 1. + baltp(imid,j,1) / re
        print *,j,preal(j)
    enddo

! read in phi from sami3

    do k = 1,ntime
        read(498) temp
        do j = 1,nny
            do i = 1,nnx
                phi(i,j,k) = temp(i,j)
            enddo
        enddo
    enddo

! calculate phi_corot

    do j = 1,nny
        phi_corot(j) = -92.e3 / 300. / preal(j) * bmag / .31
        print *,'pco',j,phi_corot(j)
    enddo

! calculate phi_nocorot

    do k = 1,ntime
        do j = 1,nny
            do i = 1,nnx
                phi_nocorot(i,j,k) = phi(i,j,k) - phi_corot(j)
            enddo
        enddo
    enddo

! write out phi_nocorot to file

    write(499) phi_nocorot

    close(497)
    close(498)
    close(499)

    END PROGRAM
