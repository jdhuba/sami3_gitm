!     *******************************************
!     *******************************************
     
!                  SAMI3_MPI-2.00

!     *******************************************
!     *******************************************

    program main

    use parameter_mod
    use variable_mod
    use namelist_mod
    use message_passing_mod
    use chemistry_mod
    use time_mod
    use atomic_mod
    use conductance_mod
    use exb_mod
    use misc_mod
    use grid_mod
    use gitm_mod


!     Some local variables

    real :: denitmp(nz,nf,nl), titmp(nz,nf,nl), vsitmp(nz,nf,nl)
    real :: deni_mnptmp(nz,nion),ti_mnptmp(nz,nion),te_mnptmp(nz)
    real :: denntmp(nz,nf,nl)
    real :: phi(nnx,nny)
    logical :: tflag,ttflag

! allocatable total matrices

! Total

    real, dimension(:,:,:,:), allocatable :: denit,dennt
    real, dimension(:,:,:,:), allocatable :: vsit, sumvsit
    real, dimension(:,:,:,:), allocatable :: tit
    real, dimension(:,:,:), allocatable   :: ut, vt, vpit
    real, dimension(:,:,:), allocatable   :: tet, tnt

!     height integrated pedersen/hall conductivities

    real, dimension(:,:,:), allocatable :: u1t, u2t, u3t, u4t
    real, dimension(:,:,:), allocatable :: vnqt, vnpt, vnphit
    real, dimension(:,:,:), allocatable :: jpt, jphit
    real, dimension(:,:,:), allocatable :: u1pt, u2st, u3ht
    real, dimension(:,:,:), allocatable :: sigmapict,sigmahict
    real, dimension(:,:,:), allocatable :: sigmapt,sigmaht

    real, dimension(:), allocatable :: dtc

! Output matrices for restart

    real :: deniout(nz,nf,nl,nion,numwork), &
            tiout(nz,nf,nl,nion,numwork), &
            vsiout(nz,nf,nl,nion,numwork)
    real :: teout(nz,nf,nl,numwork)
    real*8 :: dphi(nnx+1,nnyt)

    logical :: lprnt

!     Begin MPI stuff

    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE)
    integer :: left,right
    integer :: nerrorcode


! ************************ initializations ***********************************
! Find out how many tasks are in this partition and what my task id is.  Then
! define the number of worker tasks and the array partition size as chunksize.
! Note:  For this example, the MP_PROCS environment variable should be set
! to an odd number...to insure even distribution of the array to numtasks-1
! worker tasks.
! *****************************************************************************

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)
    write(*,*)'taskid =',taskid

    numworkers = numtasks-1

! Check to see if the number of processors selected agrees with
! the number of divisions in params

    if(taskid == 0) then
        if(numwork /= numworkers) then
            print *, ' numworkers is ',numworkers
            print *, ' numwork (in parameter_mod_mpi) is',numwork
            print *, ' in order for the code to work correctly '
            print *, ' these two numbers must be the same '
            print *, ' Either set np = numwork +1 and rerun or '
            print *, ' change numwork and recompile '
            call mpi_abort(MPI_COMM_WORLD, nerrorcode, ierr)
            call mpi_finalize(ierr)
        endif
    endif

! Determine what is left (down) and right (up)
! Here we assume that taskid=0 is the Master and does nothing but
! deal with handling the data

    print *, ' taskid ',taskid

    if(taskid == numtasks -1) then
        right = 1
    else
        right = taskid +1
    endif
    if(taskid == 1) then
        left = numtasks -1
    else
        left = taskid -1
    endif

! open daily fism file

!      if ( taskid .eq. 0 ) then
!        open(unit=402,file='fism_daily.inp ',form='unformatted')
!        read(402) fism
!        close(402)
!      endif

! open input files
! only need these on the master

    if(taskid == 0) then
        open ( unit=10, file='sami3-3.00.namelist'  )
        open ( unit=20, file='deni-init.inp'        )
        open ( unit=30, file='ichem.inp'            )
        open ( unit=50, file='phabsdt_euvac.inp'    )  ! euvac
        open ( unit=60, file='phiondt_euvac.inp'    )  ! euvac
        open ( unit=61, file='phionnt.inp'          )
        open ( unit=65, file='euvflux_euvac.inp'    )  ! euvac
        open ( unit=66, file='thetant.inp'          )
        open ( unit=67, file='zaltnt.inp'           )
    endif

    call initial

! We are out of initial now.
! We will overwrite the values of
! deni, vsi, ti, te if this is a restart (restart = true)

    if(restart) then
        if(taskid == 0) then
            print *,'doing restart'
            open ( unit=210, file='time.rst', form='formatted' )
            open ( unit=211, file='deni.rst', form='unformatted' )
            open ( unit=212, file='vsi.rst', form='unformatted' )
            open ( unit=213, file='ti.rst', form='unformatted' )
            open ( unit=214, file='te.rst', form='unformatted' )

            read(210,*) hrinit
            read(211) deniout
            read(212) vsiout
            read(213) tiout
            read(214) teout

            close (210)
            close (211)
            close (212)
            close (213)
            close (214)

            do iwrk = 1,numworkers
                do nntmp = 1,nion
                    do ktmp = 1,nl
                        do jtmp = 1,nf
                            do itmp = 1,nz
                                denitmp(itmp,jtmp,ktmp) &
                                =  deniout(itmp,jtmp,ktmp,nntmp,iwrk)
                                titmp(itmp,jtmp,ktmp) &
                                =  tiout(itmp,jtmp,ktmp,nntmp,iwrk)
                                vsitmp(itmp,jtmp,ktmp) &
                                =  vsiout(itmp,jtmp,ktmp,nntmp,iwrk)
                            enddo
                        enddo
                    enddo

                    call mpi_send(denitmp, nz*nf*nl, MPI_REAL, iwrk, 9, &
                    MPI_COMM_WORLD, ierr)
                    call mpi_send(titmp, nz*nf*nl, MPI_REAL, iwrk, 9, &
                    MPI_COMM_WORLD, ierr)
                    call mpi_send(vsitmp, nz*nf*nl, MPI_REAL, iwrk, 9, &
                    MPI_COMM_WORLD, ierr)
                enddo

                do ktmp = 1,nl
                    do jtmp = 1,nf
                        do itmp = 1,nz
                            te(itmp,jtmp,ktmp) &
                            =  teout(itmp,jtmp,ktmp,iwrk)
                        enddo
                    enddo
                enddo
                call mpi_send(te, nz*nf*nl, MPI_REAL, iwrk, 9, &
                MPI_COMM_WORLD, ierr)

            enddo
        endif

    ! Now let's get those files

        if(taskid > 0 .AND. taskid <= numworkers) then

            do nntmp = 1,nion
                call mpi_recv(denitmp, nz*nf*nl, MPI_REAL, 0, 9, &
                MPI_COMM_WORLD, status, ierr)
                call mpi_recv(titmp, nz*nf*nl, MPI_REAL, 0, 9, &
                MPI_COMM_WORLD, status, ierr)
                call mpi_recv(vsitmp, nz*nf*nl, MPI_REAL, 0, 9, &
                MPI_COMM_WORLD, status, ierr)
                do ktmp = 1,nl
                    do jtmp = 1,nf
                        do itmp = 1,nz
                            deni(itmp,jtmp,ktmp,nntmp) &
                            =  denitmp(itmp,jtmp,ktmp)
                            ti(itmp,jtmp,ktmp,nntmp) &
                            =  titmp(itmp,jtmp,ktmp)
                            vsi(itmp,jtmp,ktmp,nntmp) &
                            =  vsitmp(itmp,jtmp,ktmp)
                        enddo
                    enddo
                enddo
            enddo
            call mpi_recv(te, nz*nf*nl, MPI_REAL, 0, 9, &
            MPI_COMM_WORLD, status, ierr)

        endif

    ! tell the workers the starting time
    ! this call has to be seen by the master and workers

        call mpi_bcast(hrinit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    endif

! open output files

    if(taskid == 0) call open_u

    if(taskid == 0) then
        close (10)
        close (20)
        close (30)
        close (50)
        close (60)
        close (61)
        close (65)
        close (66)
        close (67)
        close (68)
    endif

!******************* master task *******************************************
          
    if(taskid == 0) then

    ! allocate the total matrices only on master

        allocate (denit(nz,nf,nlt,nion),dennt(nz,nf,nlt,nneut))
        allocate  (vsit(nz,nf,nlt,nion),sumvsit(nz,nf,nlt,nion))
        allocate  (tit(nz,nf,nlt,nion))
        allocate (ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt))
        allocate (tet(nz,nf,nlt),tnt(nz,nf,nlt))
        allocate (u1t(nz,nf,nlt),u2t(nz,nf,nlt),u3t(nz,nf,nlt), &
                 u4t(nz,nf,nlt))
        allocate (vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt))
        allocate (jpt(nz,nf,nlt),jphit(nz,nf,nlt))
        allocate (u1pt(nz,nf,nlt),u2st(nz,nf,nlt),u3ht(nz,nf,nlt))
        allocate (sigmapict(nz,nf,nlt),sigmahict(nz,nf,nlt))
        allocate (sigmapt(nz,nf,nlt),sigmaht(nz,nf,nlt))

        allocate (dtc(numworkers))

        hrut    = hrinit
        timemax = hrmax * sphr
        istep   = 0
        tneut   = 0.
        time    = 0.
        ntm     = 0
        ntmmax  = min(maxstep, int( hrmax / dthr ))

        print *,'max,ntmmax ',max(dt0/3600.,dthr),ntmmax

        ifintot  = numworkers
        ifintot1 = numworkers
        ifintot2 = numworkers

        tflag  = .TRUE. 
        icnt10 =  0

        do while ( tflag )

            do iwrk = 1,numworkers
                call mpi_iprobe(iwrk,10,MPI_COMM_WORLD,flagit10, &
                status,ierr)
                if (flagit10) then
                    icnt10 = icnt10 + 1
                    call mpi_recv(xxx,1,MPI_REAL,iwrk,10, &
                    MPI_COMM_WORLD,status,ierr)
                endif
                if (icnt10 == numworkers) tflag= .FALSE. 
            enddo

        ! Now wait to receive back the results from each worker task

            do  iwrk = 1, numworkers
                source = iwrk
                dest = source

                call mpi_iprobe(source, 2, &
                MPI_COMM_WORLD, flagit, status, ierr)

                if(flagit .AND. ifintot2 > 0) then

                    call mpi_recv(hipcp, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hihcm, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hipcphi, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidphig, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidpg, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidphiv, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidpv, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hipc, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hihc, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidv, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)

                    call mpi_recv(hrut, 1, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            hipcpt(j,kk)   = hipcp(j,k)
                            hihcmt(j,kk)   = hihcm(j,k)
                            hipcphit(j,kk) = hipcphi(j,k)
                            hidphigt(j,kk) = hidphig(j,k)
                            hidpgt(j,kk)   = hidpg(j,k)
                            hidphivt(j,kk) = hidphiv(j,k)
                            hidpvt(j,kk)   = hidpv(j,k)
                            hipct(j,kk)    = hipc(j,k)
                            hihct(j,kk)    = hihc(j,k)
                            hidvt(j,kk)    = hidv(j,k)
                        enddo
                    enddo


                    ifintot2 = ifintot2 - 1

                endif

                if ( ifintot2 == 0 ) then
                    ifintot2 = numworkers
                    call potpphi(phi,dphi,hrut)
                    do jwrk = 1,numworkers
                        call mpi_send(phi,nnx*nny,MPI_REAL,jwrk,3, &
                        MPI_COMM_WORLD,ierr)
                    enddo
                endif

                call mpi_iprobe(source, 0, &
                MPI_COMM_WORLD, flagit, status, ierr)

                if(flagit .AND. ifintot > 0) then
                                      
                ! This is just for outputting the data
                ! only sent as often as data dumps are requested

                    call mpi_recv(time, 1, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hrut, 1, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(istep, 1, MPI_INTEGER, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    do nntmp = 1,nion
                        call mpi_recv(denitmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(denntmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(titmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(vsitmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        do itmp = 1,nz
                            do jtmp = 1,nf
                                do ktmp = 1,nl
                                    deni(itmp,jtmp,ktmp,nntmp) &
                                      =  denitmp(itmp,jtmp,ktmp)
                                    denn(itmp,jtmp,ktmp,nntmp) &
                                      =  denntmp(itmp,jtmp,ktmp)
                                    ti(itmp,jtmp,ktmp,nntmp) &
                                      =  titmp(itmp,jtmp,ktmp)
                                    vsi(itmp,jtmp,ktmp,nntmp) &
                                      =  vsitmp(itmp,jtmp,ktmp)
                                enddo
                            enddo
                        enddo
                    enddo
                    call mpi_recv(te, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u1p, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u2s, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u3h, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u1, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u2, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u3, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u4, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmap, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmah, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmapic, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmahic, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(vnq, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(vnp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(vnphi, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(jp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(jphi, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do nn = 1,nion
                        do k = 2,nl-1
                            kk = (iwrk-1)*(nl -2) +k -1
                            if(kk == 0) kk = nlt
                            if(kk == nltp1) kk = 1
                            do j = 1,nf
                                do i = 1,nz
                                    denit(i,j,kk,nn) = deni(i,j,k,nn)
                                    dennt(i,j,kk,nn) = denn(i,j,k,nn)
                                    tit(i,j,kk,nn)   = ti(i,j,k,nn)
                                    vsit(i,j,kk,nn)  = vsi(i,j,k,nn)
                                enddo
                            enddo
                        enddo

                    ! Put the submatrices into the total matrix for restart

                        do k = 1,nl
                            do j = 1,nf
                                do i = 1,nz
                                    deniout(i,j,k,nn,iwrk) = deni(i,j,k,nn)
                                    tiout(i,j,k,nn,iwrk)   = ti(i,j,k,nn)
                                    vsiout(i,j,k,nn,iwrk)  = vsi(i,j,k,nn)
                                enddo
                            enddo
                        enddo
                    enddo  ! for nion loop

                    do k = 1,nl
                        do j = 1,nf
                            do i = 1,nz
                                teout(i,j,k,iwrk) = te(i,j,k)
                            enddo
                        enddo
                    enddo

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                tet(i,j,kk)       = te(i,j,k)
                                u1pt(i,j,kk)      = u1p(i,j,k)
                                u2st(i,j,kk)      = u2s(i,j,k)
                                u3ht(i,j,kk)      = u3h(i,j,k)
                                u1t(i,j,kk)       = u1(i,j,k)
                                u2t(i,j,kk)       = u2(i,j,k)
                                u3t(i,j,kk)       = u3(i,j,k)
                                u4t(i,j,kk)       = u4(i,j,k)
                                sigmapt(i,j,kk)   = sigmap(i,j,k)
                                sigmaht(i,j,kk)   = sigmah(i,j,k)
                                sigmapict(i,j,kk) = sigmapic(i,j,k)
                                sigmahict(i,j,kk) = sigmahic(i,j,k)
                                vnqt(i,j,kk)      = vnq(i,j,k)
                                vnpt(i,j,kk)      = vnp(i,j,k)
                                vnphit(i,j,kk)    = vnphi(i,j,k)
                                jpt(i,j,kk)       = jp(i,j,k)
                                jphit(i,j,kk)     = jphi(i,j,k)
                            enddo
                        enddo
                    enddo
                                     
                    ifintot = ifintot -1

                endif

                call mpi_iprobe(source, 1, &
                MPI_COMM_WORLD, flagit1, status, ierr)
                   
                if(flagit1 .AND. ifintot1 > 0) then
                    call mpi_recv(dtmp, 1, MPI_REAL, iwrk, 1, &
                                  MPI_COMM_WORLD, status, ierr)
                    dtc(iwrk) = dtmp
                    call mpi_recv(time, 1, MPI_REAL, iwrk, 1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(istep, 1, MPI_INTEGER, iwrk, 1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(deni_mnptmp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(ti_mnptmp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(te_mnptmp,nz,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, status, ierr)

                    if ( ifintot1 == numworkers ) then
                        do ni = nion1,nion2
                            do i = 1,nz
                                deni_mnp(i,ni) = 0.
                                ti_mnp(i,ni)   = 0.
                            enddo
                        enddo

                        do i = 1,nz
                            te_mnp(i)     = 0.
                        enddo
                    endif

                    do ni = nion1,nion2
                        do i = 1,nz
                            deni_mnp(i,ni) = deni_mnp(i,ni) + &
                            deni_mnptmp(i,ni)/numworkers
                            ti_mnp(i,ni)   = ti_mnp(i,ni) + &
                            ti_mnptmp(i,ni)/numworkers
                        enddo
                    enddo

                    do i = 1,nz
                        te_mnp(i)     = te_mnp(i) + te_mnptmp(i)/numworkers
                    enddo

                    ifintot1 = ifintot1 - 1

                endif

            enddo               ! end worker loop

        ! if we are here, we should have gathered up all the data

            if(ifintot == 0) then
                ifintot = numworkers
                ntm = ntm + 1
                call output ( hrut,ntm,istep,phi,denit,dennt,vsit, &
                              sumvsit,tit,ut,vt,vpit,tet,tnt,u1t, &
                              u2t,u3t,u4t,vnqt,vnpt,vnphit,jpt,jphit, &
                              u1pt,u2st,u3ht,sigmapict,sigmahict, &
                              sigmapt,sigmaht )
                                         
            ! write the restart files and close those files

                open ( unit=210,  file='time.rst', form='formatted' )
                open ( unit=211,  file='deni.rst', form='unformatted' )
                open ( unit=212,  file='vsi.rst', form='unformatted' )
                open ( unit=213,  file='ti.rst', form='unformatted' )
                open ( unit=214,  file='te.rst', form='unformatted' )
                open ( unit= 2322,file='dphi.rst',form='unformatted')

                write(210,*) hrut
                write(211)   deniout
                write(212)   vsiout
                write(213)   tiout
                write(214)   teout
                write(2322)  dphi

                close (210)
                close (211)
                close (212)
                close (213)
                close (214)
                close (2322)

            endif

        ! Need to fix up dt calculation

            if(ifintot1 == 0) then
                ifintot1 = numworkers
                dt = minval(dtc)
                print *,'dt = ',dt
                do  iwrk = 1,numworkers
                    call mpi_send(dt, 1, MPI_REAL, iwrk, 1, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(deni_mnp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(ti_mnp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(te_mnp,nz,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, ierr)
                enddo
            endif

        enddo       ! end while (tflag)
              
        print *, 'MASTER: All Done!'

    ! close files

        close (10)
        close (20)
        close (40)
        close (70)
        close (71)
        close (72)
        close (73)
        close (74)
        close (75)
        close (78)
        close (79)
        close (80)
        close (90)
        close (91)
        close (92)
        close (93)
        close (94)
        close (95)
        close (96)
        close (97)
        close (98)
        close (81)
        close (82)
        close (83)
        close (84)
        close (85)
        close (86)
        close (87)
        close (88)
        close (711)
        close (712)
        close (713)
        close (714)
        close (715)
        close (1712)
        close (1713)
        close (1714)
        close (1715)
        close (569)
        close (716)
        close (717)
        close (1718)
        close (811)
        close (812)
        close (813)
        close (814)
        close (815)
        close (816)
        close (817)
        close (911)
        close (912)
        close (913)
        close (914)
        close (915)
        close (916)
        close (917)
        close (384)
        close (385)
        close (386)
        close (196)
        close (197)
        close (198)
        close (201)
        close (202)
        close (491)
        close (492)
        close (493)
        close (494)
        close (495)
        close (496)
        close (497)
        close (498)

    endif

!******************* end master task ***************************************

!******************** worker task *******************************************

    if(taskid > 0) then

    ! field line loop: actual run

        hrut    = hrinit
        timemax = hrmax * sphr
        istep   = 0
        tneut   = 0.
        time    = 0.
        ttflag  = .TRUE. 
        ntm     = 0
        nfile   = 2

        ntmmax  = min(maxstep, int(hrmax / dthr))
        print *,'iwrk ',taskid,ntmmax
                 
    ! initialize neutrals
    ! neutral density, temperature, and neutral wind
    ! already done in initialization

        if (restart) then
            do nll = 1,nl
                call neutambt (hrinit,nll)
            enddo
        endif

        do while ( istep <= maxstep .AND. ttflag )

    ! parallel transport

            do nll = 1,nl
                do nfl = 1,nf
                    call zenith (hrut,nfl,nll)
                    call transprt (nfl,nll)
                enddo
            enddo

    ! Do data exchanges between guard cells

    ! buffer and send to the LEFT

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        tl1s(i,j,k) = deni(i,j,2,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        tl1s(i,j,k) = ti(i,j,2,k-nion)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    tl1s(i,j,k) = te(i,j,2)
                enddo
            enddo

            call mpi_sendrecv(tl1s, (nion+nion+1)*nz*nf, MPI_REAL, &
            left, 0, tl1r, (nion+nion+1)*nz*nf, MPI_REAL, &
            right, 0, MPI_COMM_WORLD, status, ierr)

        ! Now everybody receives

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        deni(i,j,nl,k) = tl1r(i,j,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        ti(i,j,nl,k-nion) = tl1r(i,j,k)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    te(i,j,nl) = tl1r(i,j,k)
                enddo
            enddo

        ! Buffer and send to the RIGHT

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        tr1s(i,j,k) = deni(i,j,nl-1,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        tr1s(i,j,k) = ti(i,j,nl-1,k-nion)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    tr1s(i,j,k) = te(i,j,nl-1)
                enddo
            enddo

            call mpi_sendrecv(tr1s, (nion+nion+1)*nz*nf, MPI_REAL, &
            right, 0, tr1r, (nion +nion +1)*nz*nf, MPI_REAL, &
            left, 0, MPI_COMM_WORLD, status, ierr)

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        deni(i,j,1,k) = tr1r(i,j,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        ti(i,j,1,k-nion) = tr1r(i,j,k)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    te(i,j,1) = tr1r(i,j,k)
                enddo
            enddo

        ! We are now finished exchanging guard cell data

        ! Sending hipcp and hidphig to master to calculate the
        ! potential

            call mpi_send(hipcp, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hihcm, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hipcphi, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidphig, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidpg, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidphiv, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidpv, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hipc, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hihc, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidv, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)

        ! send master current time: hrut

            call mpi_send(hrut, 1, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)

        ! now get the potential from master

            call mpi_recv(phi, nnx*nny, MPI_REAL, 0, 3, &
                          MPI_COMM_WORLD, status, ierr)

        ! perpendicular transport


            call exb  (hrut,phi)
            call neut (hrut)
            call courant

        !        print *,'out of courant, dt, taskid',dt,taskid

        !       average magnetic pole grid values (deni,Ti,Te)

            j0           = nf

            do ni = nion1,nion2
                do i = 1,nz
                    deni_mnp0 = 0.
                    ti_mnp0   = 0.
                    do k = 2,nl-1
                        if ( alts (i,j0,k) < alt_crit_avg) then
                            deni_mnp0      = deni_mnp0 + deni(i,j0,k,ni)
                            ti_mnp0        = ti_mnp0 + ti(i,j0,k,ni)
                        endif
                    enddo
                    deni_mnp(i,ni) = deni_mnp0 / float(nl-2)
                    ti_mnp(i,ni)   = ti_mnp0 / float(nl-2)
                enddo
            enddo

            do i = 1,nz
                te_mnp0 = 0.
                do k = 2,nl-1
                    if ( alts (i,j0,k) < alt_crit_avg) then
                        te_mnp0     = te_mnp0 + te(i,j0,k)
                    endif
                enddo
                te_mnp(i)   = te_mnp0 / float(nl-2)
            enddo


        ! send local dt

            call mpi_send(dt, 1, MPI_REAL, 0, 1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(time, 1, MPI_REAL, 0, 1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(istep, 1, MPI_INTEGER, 0, 1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(deni_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(ti_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(te_mnp,nz,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, ierr)

        ! get global dt

            call mpi_recv(dt, 1, MPI_REAL, 0, 1, &
                          MPI_COMM_WORLD, status, ierr)
            call mpi_recv(deni_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, status, ierr)
            call mpi_recv(ti_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, status, ierr)
            call mpi_recv(te_mnp,nz,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, status, ierr)

        ! update neutrals

        ! timing is setup up for eclipse run
        ! start at 1600 UT

        ! this is for msis

!!$            if( hrut < 16. .and. tneut >= 0.25 ) then
!!$                do nll = 2,nl-1
!!$                    call neutambt (hrut,nll)
!!$                enddo
!!$                tneut = 0.
!!$             endif

        ! this is for gitm
        ! (this is not good - running msis prior to this)

            if( tneut >= dt_gitm ) then
                call read_gitm_for_sami3(nfile)
                do nll = 1,nl
                    call neutambt (hrut,nll)
                enddo
                nfile = nfile + 1
                tneut = 0.
                if ( hrut < hrpr+hrinit ) then
                    print *,'No output yet: hr = ',hrut
                endif
            endif

        ! output data

            if ( lprnt .AND. hrut >= hrpr+hrinit) then

            ! We no longer call output from here, but send data to the MASTER
            ! The four things we want to send are  deni, ti, vsi, te

                ntm    = ntm + 1

                call mpi_send(time, 1, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(hrut, 1, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(istep, 1, MPI_INTEGER, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                do nntmp = 1,nion
                    do itmp = 1,nz
                        do jtmp = 1,nf
                            do ktmp = 1,nl
                                denitmp(itmp,jtmp,ktmp) &
                                  = deni(itmp,jtmp,ktmp,nntmp)
                                denntmp(itmp,jtmp,ktmp) &
                                  = denn(itmp,jtmp,ktmp,nntmp)
                                titmp(itmp,jtmp,ktmp) &
                                  = ti(itmp,jtmp,ktmp,nntmp)
                                vsitmp(itmp,jtmp,ktmp) &
                                  = vsi(itmp,jtmp,ktmp,nntmp)
                            enddo
                        enddo
                    enddo
                    call mpi_send(denitmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(denntmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(titmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(vsitmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                enddo

                call mpi_send(te, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(u1p, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u2s, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u3h, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u1, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u2, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u3, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u4, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(sigmap, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(sigmah, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(sigmapic, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(sigmahic, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(vnq, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(vnp, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(vnphi, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(jp, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(jphi, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                 
                dt      = dt_old
                lprnt   = .FALSE. 
                if ( taskid == 1 ) print *,'out of print dt = ',dt
                if ( ntm >= ntmmax ) ttflag = .FALSE. 

            endif

        ! time/step advancement

            istep  = istep + 1
            time   = time  + dt
            hrut   = time / sphr + hrinit
            tneut  = tneut + dt / sphr
            dts    = min(dt,mod(time,dthr*3600.))
            dt_old = dt
            if ( dts < dt .AND. .not. lprnt ) then
                dt     = max(dts,.01)
                lprnt  = .TRUE. 
            endif

        enddo    ! end time loop

    endif ! end worker task

    xxx = 1.
    call mpi_send(xxx, 1, MPI_REAL, 0, 10, &
    MPI_COMM_WORLD, ierr)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call mpi_finalize(ierr)
    print *,'done finalizing,taskid',taskid

    stop
    END PROGRAM

!******************************************
!******************************************

!            initial

!******************************************
!******************************************

    subroutine initial

    use parameter_mod
    use variable_mod
    use namelist_mod
    use message_passing_mod
    use chemistry_mod
    use time_mod
    use photo_production_mod
    use atomic_mod
    use exb_mod
    use misc_mod
    use grid_mod
 
    include 'mpif.h'

    real, dimension(:,:,:), allocatable :: altstmp,glatstmp,glonstmp
    real, dimension(:,:,:), allocatable :: altst,glatst,glonst
    real, dimension(:,:,:), allocatable :: baltst,blatst,blonst
    real, dimension(:,:,:), allocatable :: xst,yst,zst
    real, dimension(:,:,:), allocatable :: altptmp,blatptmp,blonptmp
    real, dimension(:,:,:), allocatable :: baltpt
    real, dimension(:,:,:), allocatable :: vpsnxt,vpsnyt,vpsnzt
    real, dimension(:,:,:), allocatable :: vhsnxt,vhsnyt,vhsnzt
    real, dimension(:,:,:), allocatable :: xpt,ypt,zpt
    real, dimension(:,:,:), allocatable :: bdirsxt,bdirsyt,bdirszt
    real, dimension(:,:,:), allocatable :: gsthetaxt
    real, dimension(:,:,:), allocatable :: gsthetayt
    real, dimension(:,:,:), allocatable :: gsthetazt
    real, dimension(:,:,:), allocatable :: gsphixt
    real, dimension(:,:,:), allocatable :: gsphiyt
    real, dimension(:,:,:), allocatable :: gsphizt
    real, dimension(:,:,:), allocatable :: gsrxt
    real, dimension(:,:,:), allocatable :: gsryt
    real, dimension(:,:,:), allocatable :: gsrzt
    real, dimension(:,:,:), allocatable :: xrgt
    real, dimension(:,:,:), allocatable :: xthgt
    real, dimension(:,:,:), allocatable :: xphigt

!      real fism(linesuv)

    real    :: f1026(nz,nf,nl,91),f584(nz,nf,nl,91), &
               f304 (nz,nf,nl,91),f1216(nz,nf,nl,91)
    integer :: status(MPI_STATUS_SIZE)
    real    :: zi(29),denii(29,7)
    real    :: phionr(linesuv,5)
    real    :: fluxdat(linesuv,2) ! original euvac stuff

! Some local variables

    namelist / go / maxstep,hrmax,dthr,hrpr,dt0, &
      rmin, altmin, fbar,f10p7,ap, year,day,mmass, &
      nion1,nion2,hrinit,tvn0,tvexb0,ver,veh,vw, &
      gamss, alt_crit_avg, blat_max4, blat_min,&
      snn,stn,denmin,alt_crit,cqe,plat,plon, &
      psmooth,hall,restart, &
      storm_ti,storm_tf,vexb_max, &
      lmadala,lcr,lvs,lweimer,decay_time,pcrit, &
      lhwm93,lhwm14, anu_drag0

    if ( taskid == 0 ) then
        allocate (altstmp(nz,nf,nl), glatstmp(nz,nf,nl),glonstmp(nz,nf,nl))
        allocate (altst(nz,nf,nlt),  glatst(nz,nf,nlt), glonst(nz,nf,nlt))
        allocate (baltst(nz,nf,nlt), blatst(nz,nf,nlt), blonst(nz,nf,nlt))
        allocate (xst(nz,nf,nlt),    yst(nz,nf,nlt), zst(nz,nf,nlt))
        allocate (altptmp(nzp1,nfp1,nlp1), blatptmp(nzp1,nfp1,nlp1), blonptmp(nzp1,nfp1,nlp1))
        allocate (vpsnxt(nz,nf,nlt), vpsnyt(nz,nf,nlt), vpsnzt(nz,nf,nlt))
        allocate (vhsnxt(nz,nf,nlt), vhsnyt(nz,nf,nlt), vhsnzt(nz,nf,nlt))
        allocate (xpt(nzp1,nfp1,nlt), ypt(nzp1,nfp1,nlt), zpt(nzp1,nfp1,nlt))
        allocate (bdirsxt(nz,nf,nlt), bdirsyt(nz,nf,nlt), bdirszt(nz,nf,nlt))
        allocate (gsthetaxt(nz,nf,nlt), gsthetayt(nz,nf,nlt), gsthetazt(nz,nf,nlt))
        allocate (gsphixt(nz,nf,nlt), gsphiyt(nz,nf,nlt), gsphizt(nz,nf,nlt))
        allocate (gsrxt(nz,nf,nlt), gsryt(nz,nf,nlt), gsrzt(nz,nf,nlt))
        allocate (xrgt(nz,nf,nlt), xthgt(nz,nf,nlt), xphigt(nz,nf,nlt))
        allocate (baltpt(nzp1,nfp1,nlt))
    endif

! read in parameters and initial ion density data

    if(taskid == 0) then
        read(10,go)
    endif

! send the namelist data to all the other processors

    call mpi_bcast(maxstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hrmax,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dthr,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hrpr,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dt0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(altmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fbar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(f10p7,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ap,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(year,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(day,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(mmass,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nion1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nion2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tvn0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tvexb0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ver,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(veh,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vw,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gamss,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alt_crit_avg,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blat_max4,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blat_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(snn,7,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(stn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(denmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alt_crit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cqe,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(plat,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(plon,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(psmooth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hall,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(storm_ti,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(storm_tf,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vexb_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lmadala,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lcr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lvs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lweimer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(decay_time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(pcrit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    if ( .NOT. restart) &
      call mpi_bcast(hrinit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lhwm93,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lhwm14,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(anu_drag0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    dt = dt0

    ami(pthp)  = 1.
    ami(pthep) = 4.
    ami(ptnp)  = 14.
    ami(ptop)  = 16.
    ami(ptn2p) = 28.
    ami(ptnop) = 30.
    ami(pto2p) = 32.

    amn(pth)  = 1.
    amn(pthe) = 4.
    amn(ptn)  = 14.
    amn(pto)  = 16.
    amn(ptn2) = 28.
    amn(ptno) = 30.
    amn(pto2) = 32.

    alpha0(pth)  = 0.67
    alpha0(pthe) = 0.21
    alpha0(ptn)  = 1.10
    alpha0(pto)  = 0.79
    alpha0(ptn2) = 1.76
    alpha0(ptno) = 1.74
    alpha0(pto2) = 1.59

    do i = 1,7
        aap(i) = ap
    enddo

! read in initial density data

    if(taskid == 0) then
        do i = 1,29
            read(20,102) zi(i),(denii(i,j),j=1,7)
            102 format(1x,f7.1,1p7e8.1)
        enddo
    endif

    call mpi_bcast(zi,29, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(denii,29*7, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)

! read in chemistry data
! in format statement 104 need to 'hardwire' nneut (= 7)

    if(taskid == 0) then
        do k = 1,nchem
            read(30,103) (ichem(k,j),j=1,3)
            103 format(3i3)
        enddo
    endif
    call mpi_bcast(ichem,nchem*3, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)

! generate the mesh data by everybody but the Master

    if (taskid == 0) call blonp0a
    if(taskid > 0)  call grid3_mpi

! Now wait to receive back the results from each worker task

    if(taskid == 0) then
                 
        ifintot = numworkers

        do while( ifintot > 0)

            do  iwrk = 1, numworkers
                source = iwrk
                dest = source
                            
                call mpi_iprobe(source, 0, &
                MPI_COMM_WORLD, flagit, status, ierr)
                               
                if(flagit .AND. ifintot > 0) then

                !  The three things we want to receive are  altpt blatpt blonpt

                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +  k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                baltpt(i,j,kk) = altptmp(i,j,k)
                                blatpt(i,j,kk) = blatptmp(i,j,k)
                                blonpt(i,j,kk) = blonptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are  xpt ypt zpt

                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                xpt(i,j,kk) = altptmp(i,j,k)
                                ypt(i,j,kk) = blatptmp(i,j,k)
                                zpt(i,j,kk) = blonptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! We want to receive pp

                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                ppt(i,j,kk)   = altptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! The three things we want to receive are  altst glatst glonst

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                altst(i,j,kk)  = altstmp(i,j,k)
                                glatst(i,j,kk) = glatstmp(i,j,k)
                                glonst(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! The three things we want to receive are  baltst blatst blonst

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                baltst(i,j,kk) = altstmp(i,j,k)
                                blatst(i,j,kk) = glatstmp(i,j,k)
                                blonst(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! The three things we want to receive are  xst yst zst

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                xst(i,j,kk) = altstmp(i,j,k)
                                yst(i,j,kk) = glatstmp(i,j,k)
                                zst(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! The three things we want to receive are  xrg,xthg,xphig

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                xrgt(i,j,kk)   = altstmp(i,j,k)
                                xthgt(i,j,kk)  = glatstmp(i,j,k)
                                xphigt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are vpsnx vpsny vpsnz

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                vpsnxt(i,j,kk) = altstmp(i,j,k)
                                vpsnyt(i,j,kk) = glatstmp(i,j,k)
                                vpsnzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are vhsnx vhsny vhsnz

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                vhsnxt(i,j,kk) = altstmp(i,j,k)
                                vhsnyt(i,j,kk) = glatstmp(i,j,k)
                                vhsnzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are bdirsx bdirsy bdirsz

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                bdirsxt(i,j,kk) = altstmp(i,j,k)
                                bdirsyt(i,j,kk) = glatstmp(i,j,k)
                                bdirszt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are gsthetax/y/z

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                gsthetaxt(i,j,kk) = altstmp(i,j,k)
                                gsthetayt(i,j,kk) = glatstmp(i,j,k)
                                gsthetazt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are gsphix/y/z

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                gsphixt(i,j,kk) = altstmp(i,j,k)
                                gsphiyt(i,j,kk) = glatstmp(i,j,k)
                                gsphizt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo


                !  The three things we want to receive are gsrx/y/z

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                gsrxt(i,j,kk) = altstmp(i,j,k)
                                gsryt(i,j,kk) = glatstmp(i,j,k)
                                gsrzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                    ifintot = ifintot -1
                                      
                endif

            enddo                  ! end worker loop

        enddo

    ! if we are here, we should have gathered up all the data

    ! output grid data

            open ( unit=69, file='zaltu.dat',form='unformatted' )
            open ( unit=76, file='glatu.dat',form='unformatted' )
            open ( unit=77, file='glonu.dat',form='unformatted' )
            write(69) altst
            write(76) glatst
            write(77) glonst
            close(69)
            close(76)
            close(77)

            open (144,file='glons_cg.rst',form='unformatted')
            write(144) glonst
            close(144)


            open ( unit=69, file='baltu.dat',form='unformatted' )
            open ( unit=76, file='blatu.dat',form='unformatted' )
            open ( unit=77, file='blonu.dat',form='unformatted' )
            write(69) baltst
            write(76) blatst
            write(77) blonst
            close(69)
            close(76)
            close(77)

            open (144,file='blons_cg.rst',form='unformatted')
            write(144) blonst
            close(144)

            open ( unit=69, file='xsu.dat',form='unformatted' )
            open ( unit=76, file='ysu.dat',form='unformatted' )
            open ( unit=77, file='zsu.dat',form='unformatted' )
            write(69) xst
            write(76) yst
            write(77) zst
            close(69)
            close(76)
            close(77)

            open ( unit=169, file='baltpu.dat'   ,form='unformatted' )
            open ( unit=176, file='blatpu.dat'   ,form='unformatted' )
            open ( unit=177, file='blonpu.dat'   ,form='unformatted' )
            write(169) baltpt
            write(176) blatpt
            write(177) blonpt
            close(169)
            close(176)
            close(177)

            open (144,file='blonp_cg.rst',form='unformatted')
            write(144) blonpt
            close(144)

            open ( unit=169, file='xpu.dat'   ,form='unformatted' )
            open ( unit=176, file='ypu.dat'   ,form='unformatted' )
            open ( unit=177, file='zpu.dat'   ,form='unformatted' )
            write(169) xpt
            write(176) ypt
            write(177) zpt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='xrgu.dat'   ,form='unformatted' )
            open ( unit=176, file='xthgu.dat'  ,form='unformatted' )
            open ( unit=177, file='xphigu.dat' ,form='unformatted' )
            write(169) xrgt
            write(176) xthgt
            write(177) xphigt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='vpsnxu.dat'   ,form='unformatted' )
            open ( unit=176, file='vpsnyu.dat'   ,form='unformatted' )
            open ( unit=177, file='vpsnzu.dat'   ,form='unformatted' )
            write(169) vpsnxt
            write(176) vpsnyt
            write(177) vpsnzt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='vhsnxu.dat'   ,form='unformatted' )
            open ( unit=176, file='vhsnyu.dat'   ,form='unformatted' )
            open ( unit=177, file='vhsnzu.dat'   ,form='unformatted' )
            write(169) vhsnxt
            write(176) vhsnyt
            write(177) vhsnzt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='bdirsxu.dat'   ,form='unformatted' )
            open ( unit=176, file='bdirsyu.dat'   ,form='unformatted' )
            open ( unit=177, file='bdirszu.dat'   ,form='unformatted' )
            write(169) bdirsxt
            write(176) bdirsyt
            write(177) bdirszt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='gsthetaxu.dat' ,form='unformatted' )
            open ( unit=176, file='gsthetayu.dat' ,form='unformatted' )
            open ( unit=177, file='gsthetazu.dat' ,form='unformatted' )
            write(169) gsthetaxt
            write(176) gsthetayt
            write(177) gsthetazt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='gsphixu.dat'   ,form='unformatted' )
            open ( unit=176, file='gsphiyu.dat'   ,form='unformatted' )
            open ( unit=177, file='gsphizu.dat'   ,form='unformatted' )
            write(169) gsphixt
            write(176) gsphiyt
            write(177) gsphizt
            close(169)
            close(176)
            close(177)


            open ( unit=169, file='gsrxu.dat'   ,form='unformatted' )
            open ( unit=176, file='gsryu.dat'   ,form='unformatted' )
            open ( unit=177, file='gsrzu.dat'   ,form='unformatted' )
            write(169) gsrxt
            write(176) gsryt
            write(177) gsrzt
            close(169)
            close(176)
            close(177)


    endif

    100 format (1x,1p10e16.6)

! The rest of the initialization is done also for everthing but the master

    if(taskid > 0) then

    ! MS: chicrit is the zenith angle below which the Sun is visible.
    ! For points on the surface this is just pi/2, but at higher
    ! altitudes it is bigger.

        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    coschicrit(i,j,k) = cos(pie - &
                    asin( 1./ (1. + alts(i,j,k)/re) ))
                enddo
            enddo
        enddo


    ! put deni on mesh via linear interpolation
    ! and put on lower limit

    ! initialize all ions

        j0 = 1
        do n = 1,nion
            do k = 1,nl
                do j = 1,nf
                    do i = 1,nz
                        jj = 1
                        do while (  alts(i,j,k) >= zi(jj) .AND. jj <= 28 )
                            j0 = jj
                            jj = jj + 1
                        enddo
                        if ( n == 1 ) nn = pthp
                        if ( n == 2 ) nn = pthep
                        if ( n == 3 ) nn = ptnp
                        if ( n == 4 ) nn = ptop
                        if ( n == 5 ) nn = ptn2p
                        if ( n == 6 ) nn = ptnop
                        if ( n == 7 ) nn = pto2p
                        slope   = ( denii(j0+1,n) - denii(j0,n) ) &
                                / ( zi   (j0+1)   - zi   (j0) )
                        deni(i,j,k,nn) = denii(j0,n) + &
                                       ( alts(i,j,k) - zi(j0) ) * slope
                        deni(i,j,k,nn) = amax1 ( 3. * deni(i,j,k,nn) , denmin )
                        if ( alts(i,j,k) > zi(29) ) then
                            if ( n == 1 )  then
                                nn = pthp
                                deni(i,j,k,nn) = &
                                   amax1(denii(29,n)*zi(29)/alts(i,j,k),denmin)
                            else
                                deni(i,j,k,nn) = denmin
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

    !     initialize helium density = 10% hydrogen density

        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    deni(i,j,k,pthep) = 0.1 * deni(i,j,k,pthp)
                enddo
            enddo
        enddo

    ! print *,'done initializing deni',taskid

    ! initialize neutrals
    ! neutral density, temperature, and neutral wind

        if ( .NOT. restart ) then
            call read_gitm_for_sami3(1)
            do nll = 1,nl
                call neutambt (hrinit,nll)
            enddo
        endif

    ! electron and ion temperature initialization

        do k = nion1,nion2
            do n = 1,nl
                do j = 1,nf
                    do i = 1,nz
                        ti(i,j,n,k)    = tni(i,j,n)
                    enddo
                enddo
            enddo
        enddo

        do n = 1,nl
            do j = 1,nf
                do i = 1,nz
                    te(i,j,n)      = tni(i,j,n)
                enddo
            enddo
        enddo

    !       average magnetic pole grid values (deni,Ti,Te)

        j0  = nf

        do ni = nion1,nion2
            do i = 1,nz
                deni_mnp0 = 0.
                ti_mnp0   = 0.
                do k = 2,nl-1
                    deni_mnp0      = deni_mnp0 + deni(i,j0,k,ni)
                    ti_mnp0        = ti_mnp0 + ti(i,j0,k,ni)
                enddo
                deni_mnp(i,ni) = deni_mnp0 / float(nl-2)
                ti_mnp(i,ni)   = ti_mnp0 / float(nl-2)
            enddo
        enddo

        do i = 1,nz
            te_mnp0 = 0.
            do k = 2,nl-1
                te_mnp0     = te_mnp0 + te(i,j0,k)
            enddo
            te_mnp(i)   = te_mnp0 / float(nl-2)
        enddo



    ! initialize ion velocity to zero

        do nn = nion1,nion2
            do k = 1,nl
                do j = 1,nf
                    do i = 1,nz
                        vsi(i,j,k,nn)     = 0.
                        sumvsi(i,j,k,nn)  = 0.
                    enddo
                enddo
            enddo
        enddo

    endif

! endif for taskid > 0 initialization

! read in photoabsorption rates

    if(taskid == 0) then
        do i = 1,linesuv
            read (50,105) (sigabsdt(i,j), j=1,3)
            105 format (3f7.2)
        enddo
    endif

    if(taskid == 0) then
        do j = 1,3
            do i = 1,linesuv
                sigabsdt(i,j) = tm18 * sigabsdt(i,j)
            enddo
        enddo
    endif
    call mpi_bcast(sigabsdt,linesuv*3, &
    MPI_REAL,0,MPI_COMM_WORLD,ierr)

! initialize photoionization rates to zero

    do j = 1,nneut
        do i = 1,linesuv
            sigidt(i,j)  = 0.
        enddo
        do i = 1,linesnt
            sigint(i,j)  = 0.
        enddo
    enddo

! read in daytime photoionization line data
! (only n, o, he, n_2, o_2)

    if(taskid == 0) then
        do i = 1,linesuv
            read(60,106) (phionr(i,j), j=1,5)
            sigidt(i,ptn ) = phionr(i,1)
            sigidt(i,pto ) = phionr(i,2)
        ! can increase He+ photoproduction rate
        ! bailey and sellek used 2.5
        ! JK used 1.5
            sigidt(i,pthe) = phionr(i,3)
            sigidt(i,ptn2) = phionr(i,4)
            sigidt(i,pto2) = phionr(i,5)
        enddo
        106 format(5f7.2)
                 
        do j = 1,nion
            do i = 1,linesuv
                sigidt(i,j) = tm18 * sigidt(i,j)
            enddo
        enddo
                 
    ! read in nighttime photoionization line data
    ! (only o, n_2, n0, o_2)

        do i = 1,linesnt
            read(61,106) (phionr(i,j), j=1,4)
            sigint(i,pto ) = phionr(i,1)
            sigint(i,ptn2) = phionr(i,2)
            sigint(i,ptno) = phionr(i,3)
            sigint(i,pto2) = phionr(i,4)
        enddo
                 
        do j = 1,nion
            do i = 1,linesnt
                sigint(i,j) = tm18 * sigint(i,j)
            enddo
        enddo
    endif
    call mpi_bcast(sigidt,linesuv*7, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(sigint,linesnt*7, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)

!     below is altered from original
!     now just for flux spectra from harry warren

!      if(taskid .eq. 0) then
!         do i = 1,linesuv
!           flux(i) = fism(i)
!         enddo
!       endif

     
! read in f74113, ai data and set euv flux
! (from richards et al., jgr 99, 8981, 1994)
     
    p  = 0.5 * ( f10p7 + fbar )
      
    if(taskid == 0) then
        do i = 1,linesuv
            read (65,107) (fluxdat(i,j),j=1,2)
            f74   = fluxdat(i,1)
        !             if ( i .eq. 1 ) f74 = 4.4 * f74
            ai    = fluxdat(i,2)
            flux(i) = f74 * ( 1. + ai * ( p - 80.) ) * 1.e9
        enddo
    endif
    107 format (f6.3,1pe11.4)



    call mpi_bcast(flux,linesuv, &
    MPI_REAL,0,MPI_COMM_WORLD,ierr)


! read in angles for nighttime deposition fluxes

    if(taskid == 0) then
        do i = 1,linesnt
            read(66,108) (thetant(i,j), j=1,4)
        enddo
    endif
    108 format (4f7.1)
    call mpi_bcast(thetant,linesnt*4, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)


! read in min/max altitude for nighttime deposition fluxes
!   zaltnt(i,1): zmin(i)
!   zaltnt(i,2): zmax(i)

    if(taskid == 0) then
        do i = 1,linesnt
            read(67,108) (zaltnt(i,j), j=1,2)
        enddo
    endif

    call mpi_bcast(zaltnt,linesnt*2, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)

! Do this for everything but the master

    if (taskid > 0) then

    ! call nighttime euv flux subroutines
    ! (lyman beta 1026, he i 584, he ii 304, lyman alpha 1216)

        do k = 1,nl
            do j = 1,nf
                call sf1026 ( f1026,1,j,k )
                call sf584  ( f584 ,2,j,k )
                call sf304  ( f304 ,3,j,k )
                call sf1216 ( f1216,4,j,k )
                do l = 1,91
                    do i = 1,nz
                        fluxnt(i,j,k,l,1) = f1026(i,j,k,l)
                        fluxnt(i,j,k,l,2) = f584 (i,j,k,l)
                        fluxnt(i,j,k,l,3) = f304 (i,j,k,l)
                        fluxnt(i,j,k,l,4) = f1216(i,j,k,l)
                    enddo
                enddo
            enddo
        enddo

    ! initialize e x b drift to 0

        do k = 1,nl
            do j = 1,nf
                do i = 1,nzp1
                    vexbs(i,j,k) = 0.
                enddo
            enddo
        enddo

        do k = 1,nl
            do j = 1,nfp1
                do i = 1,nz
                    vexbp(i,j,k) = 0.
                enddo
            enddo
        enddo

        do k = 1,nlp1
            do j = 1,nf
                do i = 1,nz
                    vexbh(i,j,k) = 0.
                enddo
            enddo
        enddo

    ! intialize diagnostic variables to 0

        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    u1p(i,j,k) = 0.
                    u2s(i,j,k) = 0.
                    u3h(i,j,k) = 0.
                    u1(i,j,k) = 0.
                    u2(i,j,k) = 0.
                    u3(i,j,k) = 0.
                    u4(i,j,k) = 0.
                    u5(i,j,k) = 0.
                enddo
            enddo
        enddo

!        do k = 1,nion
!            do n = 1,nl
!                do j = 1,nf
!                    do i = 1,nz
!                        t1(i,j,n,k) = 0.
!                        t2(i,j,n,k) = 0.
!                        t3(i,j,n,k) = 0.
!                    enddo
!                enddo
!            enddo
!        enddo
    endif

    if ( taskid == 0 ) then
        deallocate (altstmp,glatstmp,glonstmp)
        deallocate (altst,glatst,glonst)
        deallocate (baltst,blatst,blonst)
        deallocate (xst,yst,zst)
        deallocate (altptmp,blatptmp,blonptmp)
        deallocate (baltpt)
        deallocate (vpsnxt,vpsnyt,vpsnzt)
        deallocate (vhsnxt,vhsnyt,vhsnzt)
        deallocate (xpt,ypt,zpt)
        deallocate (bdirsxt,bdirsyt,bdirszt)
        deallocate (gsthetaxt,gsthetayt,gsthetazt)
        deallocate (gsphixt,gsphiyt,gsphizt)
        deallocate (gsrxt,gsryt,gsrzt)
        deallocate (xrgt,xthgt,xphigt)
    endif

    print *,' finished initialization taskid = ',taskid

    return
    end subroutine initial


!******************************************
!******************************************

!            transprt

!******************************************
!******************************************

    subroutine transprt (nfl,nll)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use chemistry_mod
    use time_mod
    use atomic_mod
    use conductance_mod
    use grid_mod

    real :: prod(nz,nion),loss(nz,nion),lossr, &
    phprodr(nz,nion),chrate(nz,nchem), &
    chloss(nz,nion),chprod(nz,nion),relossr(nz,nion)
    real :: deni_old(nz,nion),te_old(nz),ti_old(nz,nion),vsi_old(nz,nion)
    real :: tvn(nz,nl)
    real :: nuin(nz,nion,nneut), &
    nuij(nz,nion,nion),sumnuj(nz,nion)
    real :: vsin(nz,nion),vsidn(nz,nion),denin(nz,nion),cs(nz,nion)
    real :: ten(nz),tin(nz,nion)

! calculation of production and loss
!   phprodr: photo production rates
!   chrate:  chemical rates (ichem)
!   chloss:  chemical loss term
!   chprod:  chemical production term
!   relossr: recombination loss rates

! initialize tvn and gs and cfs (centrifugal force)

    do i = 1,nz
        tvn(i,nll) = 0.
        gs(i,nll)  = 0.
        cfs(i,nll)  = 0.
    enddo

    do i = 1,nz
        ne(i,nfl,nll)   = 1.
        te_old(i)       = te(i,nfl,nll)
        do j = nion1,nion2
            deni_old(i,j) = deni(i,nfl,nll,j)
            ne(i,nfl,nll) = ne(i,nfl,nll) + deni(i,nfl,nll,j)
            ti_old(i,j)   = ti(i,nfl,nll,j)
            vsi_old(i,j)  = vsi(i,nfl,nll,j)
        enddo
    enddo

    call photprod ( phprodr,nfl,nll               ) ! calculates phprodr
    call chemrate ( chrate,nfl,nll                ) ! calculates chrate
    call chempl   ( chrate,chloss,chprod,nfl,nll  ) ! calcualtes chloss,chprod
    call recorate ( relossr,nfl,nll               ) ! calculates relossr

    do i = 1,nz
        do j = nion1,nion2
            prod  (i,j) =  phprodr(i,j) * denn(i,nfl,nll,j) &
            + chprod(i,j)
            lossr       =  relossr(i,j) * deni(i,nfl,nll,j) * &
            ne(i,nfl,nll) &
            + chloss(i,j)
            loss (i,j)  =  lossr / deni(i,nfl,nll,j)
        enddo

!  add loss for NO+ below 90 km
!  NEED TO IMPROVE THIS

    if ( alts(i,nfl,nll) .lt. 90. ) then
!      prod(i,ptnop) = 0.
      loss(i,ptnop) = loss(i,pto2p)
      loss(i,pthp)  = loss(i,pto2p)
    endif

    !     loss term for hydrogen and helium

        if ( alts(i,nfl,nll) > pcrit*re ) then
            loss(i,pthp)  = loss(i,pthp)  + 1./decay_time
            loss(i,pthep) = loss(i,pthep) + 1./decay_time
        ! add loss term for O+ (JK)
        !          loss(i,ptop) = loss(i,ptop) + 1./decay_time
        endif

        gs(i,nll)   =  gzero * xrg(i,nfl,nll) &
        * ( re / (re + alts(i,nfl,nll)) ) ** 2


    !        gs(i,nll)   = -gzero
    !     .                 * ( re / (re + alts(i,nfl,nll)) ) ** 2
    !     .                 * ( gsrx(i,nfl,nll)*bdirsx(i,nfl,nll) +
    !     .                     gsry(i,nfl,nll)*bdirsy(i,nfl,nll) +
    !     .                     gsrz(i,nfl,nll)*bdirsz(i,nfl,nll)  )


    ! K     centrifugal force (see notes 2012/01/04)

        fzero = 3.369
        clat = cos(pie*glats(i,nfl,nll)/180.0)
        slat = sin(pie*glats(i,nfl,nll)/180.0)
        cfs(i,nll)   =  -fzero * &
        (clat*xrg(i,nfl,nll) + slat*xthg(i,nfl,nll)) &
        * (re + alts(i,nfl,nll)) * clat / re


!       note: sign of gp expicitly accounted for
!       in derivation, i.e., g = -gp phat
!       so gp is positive here (JH)

        gp(i,nfl,nll)   = gzero &
        * ( re / (re + alts(i,nfl,nll)) ) ** 2 &
        * ( gsrx(i,nfl,nll)*vpsnx(i,nfl,nll) + &
        gsry(i,nfl,nll)*vpsny(i,nfl,nll) + &
        gsrz(i,nfl,nll)*vpsnz(i,nfl,nll)  )

        vnq(i,nfl,nll) = v(i,nfl,nll) * &
        ( gsthetax(i,nfl,nll) * bdirsx(i,nfl,nll) + &
        gsthetay(i,nfl,nll) * bdirsy(i,nfl,nll) + &
        gsthetaz(i,nfl,nll) * bdirsz(i,nfl,nll)   ) + &
        u(i,nfl,nll) * &
        ( gsphix(i,nfl,nll) * bdirsx(i,nfl,nll) + &
        gsphiy(i,nfl,nll) * bdirsy(i,nfl,nll) + &
        gsphiz(i,nfl,nll) * bdirsz(i,nfl,nll)   )   + &
        w(i,nfl,nll) * &
        ( gsrx(i,nfl,nll) * bdirsx(i,nfl,nll) + &
        gsry(i,nfl,nll) * bdirsy(i,nfl,nll) + &
        gsrz(i,nfl,nll) * bdirsz(i,nfl,nll)   )

        vnp(i,nfl,nll) = v(i,nfl,nll) * &
        ( gsthetax(i,nfl,nll) * vpsnx(i,nfl,nll) + &
        gsthetay(i,nfl,nll) * vpsny(i,nfl,nll) + &
        gsthetaz(i,nfl,nll) * vpsnz(i,nfl,nll)   ) + &
        u(i,nfl,nll) * &
        ( gsphix(i,nfl,nll) * vpsnx(i,nfl,nll) + &
        gsphiy(i,nfl,nll) * vpsny(i,nfl,nll) + &
        gsphiz(i,nfl,nll) * vpsnz(i,nfl,nll)   )   + &
        w(i,nfl,nll) * &
        ( gsrx(i,nfl,nll) * vpsnx(i,nfl,nll) + &
        gsry(i,nfl,nll) * vpsny(i,nfl,nll) + &
        gsrz(i,nfl,nll) * vpsnz(i,nfl,nll)   )

        vnphi(i,nfl,nll) = v(i,nfl,nll) * &
        ( gsthetax(i,nfl,nll) * vhsnx(i,nfl,nll) + &
        gsthetay(i,nfl,nll) * vhsny(i,nfl,nll) + &
        gsthetaz(i,nfl,nll) * vhsnz(i,nfl,nll)   ) + &
        u(i,nfl,nll) * &
        ( gsphix(i,nfl,nll) * vhsnx(i,nfl,nll) + &
        gsphiy(i,nfl,nll) * vhsny(i,nfl,nll) + &
        gsphiz(i,nfl,nll) * vhsnz(i,nfl,nll)   )   + &
        w(i,nfl,nll) * &
        ( gsrx(i,nfl,nll) * vhsnx(i,nfl,nll) + &
        gsry(i,nfl,nll) * vhsny(i,nfl,nll) + &
        gsrz(i,nfl,nll) * vhsnz(i,nfl,nll)   )

        tvn(i,nll)    = vnq(i,nfl,nll)
                  
    enddo

    call update ( tvn,nuin,sumnuj,nuij,nfl,nll )

    do i = 1,nz
        do nni = nion1,nion2
            sumvsi(i,nfl,nll,nni) = 0.
            do nj = nion1,nion2
                sumvsi(i,nfl,nll,nni) =   sumvsi(i,nfl,nll,nni) + &
                nuij(i,nni,nj)*vsi(i,nfl,nll,nj)
            enddo
        enddo
    enddo

! define new arrays for velocity and density

    do ni = nion1,nion2
        do i = 1,nz
            vsin (i,ni) = vsi(i,nfl,nll,ni)
            vsidn(i,ni) = vsid(i,nfl,nll,ni)
            denin(i,ni) = deni(i,nfl,nll,ni)
        enddo
    enddo

! define sound velocity used in vsisolv

    do ni = nion1,nion2
        do i = 1,nz
            cfac     = 1.6667 * 8.6174e-5 * te(i,nfl,nll) / ami(ni)
            cs(i,ni) = 9.79e5 * sqrt(cfac)
        enddo
    enddo

! update variables

    do ni = nion1,nion2

        call vsisolv ( vsin(1,ni),vsidn(1,ni),vsi_old(1,ni), &
        sumnuj(1,ni),nfl,nll,cs(1,ni) )

    ! compensating filter

        call smoothz ( vsin(1,ni), 1 )

    ! put stuff back into velocity array

        do i = 1,nz
            vsi(i,nfl,nll,ni)  = vsin(i,ni)
            vsid(i,nfl,nll,ni) = vsidn(i,ni)
        enddo

        call densolv2 ( ni,denin(1,ni), &
        prod(1,ni),loss(1,ni),deni_old(1,ni),nfl,nll )

    ! put stuff back into density array

        do i = 1,nz
            deni(i,nfl,nll,ni) = denin(i,ni)
        enddo

    ! put floor on density

        do i = 1,nz
            deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), denmin )
        ! below commented out (JK)
            if ( alts(i,nfl,nll) > pcrit*re .AND. ni == pthp ) &
            deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), .1 )
            if ( alts(i,nfl,nll) > pcrit*re .AND. ni == pthep ) &
            deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), .01 )

        enddo

    enddo


! define new arrays for temperature

    do ni = nion1,nion2
        do i = 1,nz
            tin(i,ni)  = ti(i,nfl,nll,ni)
        enddo
    enddo

    do i = 1,nz
        ten(i)  = te(i,nfl,nll)
    enddo

! temperatures (with floors and warnings)

    tn0 = 200. ! floor on temperature

    call etemp  (ten,te_old,phprodr,nfl,nll)
    do i = 1,nz
        te(i,nfl,nll)  = amax1(tn(i,nfl,nll),ten(i))
    !        te(i,nfl,nll)  = amax1(tn0,ten(i))
        te(i,nfl,nll)  = amin1(te(i,nfl,nll),1.e4)
        if ( te(i,nfl,nll) < 0 ) then
            print *,' T(e) negative: i,nfl,nll taskid',i,nfl,nll,taskid
            stop
        endif
    enddo



    call htemp  (tin(1,pthp) ,ti_old(1,pthp) ,tvn,nuin,nfl,nll)
    do i = 1,nz
        ti(i,nfl,nll,pthp)  = amax1(tn(i,nfl,nll),tin(i,pthp))
    !        ti(i,nfl,nll,pthp)  = amax1(tn0,tin(i,pthp))
        ti(i,nfl,nll,pthp)  = amin1(ti(i,nfl,nll,pthp),1.e4)
        if ( ti(i,nfl,nll,pthp) < 0 ) then
            print *,' T(H) negative: i,nfl,nll',i,nfl,nll
            stop
        endif
    enddo

    call hetemp (tin(1,pthep),ti_old(1,pthep),tvn,nuin,nfl,nll)
    do i = 1,nz
        ti(i,nfl,nll,pthep)  = amax1(tn(i,nfl,nll),tin(i,pthep))
    !        ti(i,nfl,nll,pthep)  = amax1(tn0,tin(i,pthep))
        ti(i,nfl,nll,pthep)  = amin1(ti(i,nfl,nll,pthep),1.e4)
        if ( ti(i,nfl,nll,pthep) < 0 ) then
            print *,' T(He) negative: i,nfl,nll',i,nfl,nll
            stop
        endif
    enddo

    call otemp  (tin(1,ptop) ,ti_old(1,ptop) ,tvn,nuin,nfl,nll)
    do i = 1,nz
        ti(i,nfl,nll,ptop)  = amax1(tn(i,nfl,nll),tin(i,ptop))
    !        ti(i,nfl,nll,ptop)  = amax1(tn0,tin(i,ptop))
        ti(i,nfl,nll,ptop)  = amin1(ti(i,nfl,nll,ptop),1.e4)
        if ( ti(i,nfl,nll,ptop) < 0 ) then
            print *,' T(O) negative: i,nfl,nll',i,nfl,nll
            stop
        endif
    enddo

    do i = 1,nz
        ti(i,nfl,nll,ptnp )    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptn2p)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptnop)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,pto2p)    = ti(i,nfl,nll,ptop)
    enddo

    return
    end subroutine transprt




!******************************************
!******************************************

!            update

!******************************************
!******************************************

    subroutine update ( tvn,nuin,sumnuj,nuij,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use conductance_mod
    use grid_mod

    real :: nuin(nz,nion,nneut),nuij(nz,nion,nion)
    real :: nuint(nz,nion)
    real :: sumnuj(nz,nion),nufacij,nufacin
    real :: tvn(nz,nl)
    real :: k0,mi
    real :: nuen(nz,nf,nl)

! ion-neutral collision frequency

! initialize everything to 0

    do nn = 1,nneut
        do ni = nion1,nion2
            do iz = 1,nz
                nuin (iz,ni,nn) = 0.
                nuint(iz,ni)    = 0.
            enddo
        enddo
    enddo

! collision frequencies/factors

! hydrogen (H)

    ni = pthp
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == pto ) then
                teff    = ti(i,nfl,nll,ni)
                fac     = ( 1.00 - .047 * alog10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn)  = 6.61e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
                amimn   = amn(nn) / ( ami(ni) + amn(nn) )
                nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
                nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! helium (He)

    ni = pthep
    do nn = 1,nneut
        do i = 1,nz
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! nitrogen (N)

    ni = ptnp
    do nn = 1,nneut
        do i = 1,nz
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! oxygen (O)

    ni = ptop
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == pto ) then
                teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
                fac     = ( 1.04 - .067 * alog10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn)  = 4.45e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
                amimn   = amn(nn) / ( ami(ni) + amn(nn) )
                nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
                nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! nitrogen 2(N2)

    ni = ptn2p
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == ptn2 ) then
                teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
                fac     = ( 1.00 - .069 * alog10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn) = 5.14e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
                amimn   = amn(nn) / ( ami(ni) + amn(nn) )
                nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
                nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! nitrous oxide (N0)

    ni = ptnop
    do nn = 1,nneut
        do i = 1,nz
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! oxygen 2(O2)

    ni = pto2p
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == pto2 ) then
                teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
                fac     = ( 1.00 - .073 * alog10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn) = 2.59e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
                amimn   = amn(nn) / ( ami(ni) + amn(nn) )
                nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
                nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! ion-ion collision frequency

    do ni = nion1,nion2
        do i = 1,nz
            do nj = nion1,nion2
                if(ni /= nj) then
                    alame1  = ( ami(ni) + ami(nj) ) * evtok / &
                    ( ami(ni)*ti(i,nfl,nll,nj) + &
                    ami(nj)*ti(i,nfl,nll,ni) )
                    alame2  = deni(i,nfl,nll,ni) * evtok / ti(i,nfl,nll,ni) + &
                    deni(i,nfl,nll,nj) * evtok / ti(i,nfl,nll,nj)
                    if ( alame2 < 0 ) then
                        print *,'ni,i,nj,nfl,nll,tii,tij,alame1,alame2,nii,nij', &
                        ni,i,nj,nfl,nll,ti(i,nfl,nll,ni),ti(i,nfl,nll,nj), &
                        alame1,alame2, &
                        deni(i,nfl,nll,ni),deni(i,nfl,nll,nj)
                        stop
                    endif
                    alame   = alame1 * sqrt(alame2)
                    alam    = 23. - alog(alame)
                    amufac  = (ami(nj)/ami(ni))/(ami(ni) +ami(nj))
                    nufacij = 9.2e-2*alam*sqrt(amufac)
                    nuij(i,ni,nj) =  nufacij * deni(i,nfl,nll,nj) &
                    / sqrt( ti(i,nfl,nll,ni)**3 )
                else
                    nuij(i,ni,nj) = 0.
                endif
            enddo
        enddo
    enddo

!     add this for restart (sarah mcdonald)
!     it's used in term 1 below

    do i = 1,nz
        do nni = nion1,nion2
            sumvsi(i,nfl,nll,nni) = 0.
            do nj = nion1,nion2
                sumvsi(i,nfl,nll,nni) =   sumvsi(i,nfl,nll,nni) + &
                nuij(i,nni,nj)*vsi(i,nfl,nll,nj)
            enddo
        enddo
    enddo

! sumnuj: sum of ion-ion coll freq and nuin

    do ni = nion1,nion2
        do i = 1,nz
            sumnuj(i,ni) = 0.
            do nj = nion1,nion2
                sumnuj(i,ni) = sumnuj(i,ni) + nuij(i,ni,nj)
            enddo
            sumnuj(i,ni) = sumnuj(i,ni) + nuint(i,ni)
        enddo
    enddo

! update ne

    do i = 1,nz
        ne(i,nfl,nll) = 1.
        do ni = nion1,nion2
            ne(i,nfl,nll) = ne(i,nfl,nll) + deni(i,nfl,nll,ni)
        enddo
    enddo

! get a new value for vsid

    do i = 2,nz-1
        do ni = nion1,nion2
            mi    = amu * ami(ni)
            k0    = bolt / mi
            term1 = nuint(i,ni) * tvn(i,nll) + &
            sumvsi(i,nfl,nll,ni) + gs(i,nll) + cfs(i,nll)
            pip   = 0.5 * (   deni(i+1,nfl,nll,ni) * ti(i+1,nfl,nll,ni) &
            + deni(i,nfl,nll,ni)   * ti(i,nfl,nll,ni)   )
            pim   = 0.5 * (   deni(i,nfl,nll,ni)   * ti(i,nfl,nll,ni) &
            + deni(i-1,nfl,nll,ni) * ti(i-1,nfl,nll,ni) )
            denid = &
            (        deni(i-1,nfl,nll,ni) &
            + 4. * deni(i,nfl,nll,ni) &
            +      deni(i+1,nfl,nll,ni)  ) / 6.
            term2 =  - bms(i,nfl,nll) * k0 /  denid &
            * ( pip - pim ) / d22s(i,nfl,nll)
            pep   = 0.5 * (   ne(i+1,nfl,nll) * te(i+1,nfl,nll) &
            + ne(i,nfl,nll)   * te(i,nfl,nll)   )
            pem   = 0.5 * (   ne(i,nfl,nll)   * te(i,nfl,nll) &
            + ne(i-1,nfl,nll) * te(i-1,nfl,nll) )
            dened = &
            ( ne(i-1,nfl,nll) + 4. * ne(i,nfl,nll) + ne(i+1,nfl,nll) ) / 6.
            term3 =  - bms(i,nfl,nll) * k0 /  dened &
            * ( pep - pem ) / d22s(i,nfl,nll)

            vsid(i,nfl,nll,ni)  =  term1 + term2 + term3

        enddo
    enddo

! fix up end points for vsid

    do ni = nion1,nion2
        vsid (1,nfl,nll,ni)    = vsid (2,nfl,nll,ni)
        vsid (nz,nfl,nll,ni)   = vsid (nz-1,nfl,nll,ni)
    enddo

! calculate the electron-neutral collision frequency
! nuen = 5.4e-10*n_n*T_e^1/2 (kelley, the earth's ionosphere, p. 462)

    do i = 1,nz
        nuen(i,nfl,nll) = 0
        do nn = 1,nneut
            nuen(i,nfl,nll) = nuen(i,nfl,nll) + 5.4e-10 * &
            denn(i,nfl,nll,nn) * sqrt(te(i,nfl,nll))
        enddo
    enddo


! calculate pedersen and hall conductivities

    do i = 1,nz
        dene    = ne(i,nfl,nll)
        oce     = 1.76e7 * bmag * bms(i,nfl,nll)
        sige    = dene * charge * sol / ( bmag * bms(i,nfl,nll) )
        cole    = nuen(i,nfl,nll) / oce
        denome  = 1. + cole * cole
        sigpe   = sige * cole / denome
        sighe   = sige * cole * cole / denome
        sigpi   = 0.
        sighi   = 0.
        sighic  = 0.
        sigpic  = 0.
        do ni = nion1,nion2
            oci    = 9580. * bmag * bms(i,nfl,nll) / ami(ni)
            sigi   = deni(i,nfl,nll,ni) * charge * sol / &
            ( bmag * bms(i,nfl,nll) )
            coli   = nuint(i,ni) / oci
            denomi = 1. + coli * coli
            sigpi  = sigpi  + sigi * coli / denomi
            sigpic = sigpic + sigi * coli / denomi / oci
            sighi  = sighi  + sigi * coli * coli / denomi
            sighic = sighic + sigi / denomi / oci
        enddo
        sigmap(i,nfl,nll)   = sigpi + sigpe
        sigmah(i,nfl,nll)   = sighi - sighe
        sigmapic(i,nfl,nll) = sigpic
        sigmahic(i,nfl,nll) = sighic
        if (alts(i,nfl,nll) >= 1.e4) then
            sigmap(i,nfl,nll)   = 0.00
            sigmah(i,nfl,nll)   = 0.00
            sigmapic(i,nfl,nll) = 0.00
            sigmahic(i,nfl,nll) = 0.00
        endif
    enddo

    if ( .NOT. hall ) then
        do i=1,nz
            sigmah(i,nfl,nll)   = 0.
            sigmahic(i,nfl,nll) = 0.
        enddo
    endif

    hipcp(nfl,nll)    = 0.
    hipcphi(nfl,nll)  = 0.
    hihcm(nfl,nll)    = 0.
    hidpv(nfl,nll)    = 0.
    hidphiv(nfl,nll)  = 0.
    hidpg(nfl,nll)    = 0.
    hidphig(nfl,nll)  = 0.

    hipc(nfl,nll)     = 0.
    hihc(nfl,nll)     = 0.
    hidv(nfl,nll)     = 0.
         
    do i = 1,nz
        ang     = .5 * pie - blats(i,nfl,nll) * pie / 180.
        bang    = blats(nz,nfl,nll) * pie / 180.
        del     = sqrt ( 1. + 3. * cos(ang) * cos(ang) )
        b       = bmag * bms(i,nfl,nll)

        hipcp(nfl,nll) = hipcp(nfl,nll) + &
        sigmap(i,nfl,nll) * del / bms(i,nfl,nll) * &
        dels(i,nfl,nll) 

        hipcphi(nfl,nll) = hipcphi(nfl,nll) + &
        sigmap(i,nfl,nll) / del / bms(i,nfl,nll) * &
        dels(i,nfl,nll) 

        hihcm(nfl,nll) = hihcm(nfl,nll) + &
        sigmah(i,nfl,nll) / bms(i,nfl,nll) * &
        dels(i,nfl,nll)

        fdpv = (bmag/sol) * ( sigmap(i,nfl,nll) * vnphi(i,nfl,nll) + &
        sigmah(i,nfl,nll) * vnp(i,nfl,nll)     )

        hidpv(nfl,nll) = hidpv(nfl,nll) + &
        brs(i,nfl,nll) * 1.e5  * sin(ang) * &
        fdpv * dels(i,nfl,nll)

        fdpg = (bmag/sol) * sigmapic(i,nfl,nll) * gp(i,nfl,nll)

        hidpg(nfl,nll) = hidpg(nfl,nll) + &
        brs(i,nfl,nll) * 1.e5  * sin(ang) * &
        fdpg * dels(i,nfl,nll)

        fdphiv = (bmag/sol) * ( -sigmap(i,nfl,nll) * vnp(i,nfl,nll) + &
        sigmah(i,nfl,nll) * vnphi(i,nfl,nll)  )

        hidphiv(nfl,nll) = hidphiv(nfl,nll) + &
        re * 1.e5 * ( sin(ang) ** 3 ) / del * &
        fdphiv * dels(i,nfl,nll)

        fdphig = (bmag/sol) * sigmahic(i,nfl,nll) * gp(i,nfl,nll)

        hidphig(nfl,nll) = hidphig(nfl,nll) + &
        re * 1.e5 * ( sin(ang) ** 3 ) / del * &
        fdphig * dels(i,nfl,nll)

    !       integrated quantities for current

        hipc(nfl,nll) = hipc(nfl,nll) + &
        sigmap(i,nfl,nll) * del / re / sin(ang) ** 3 * &
        dels(i,nfl,nll) / 1.e5
        hihc(nfl,nll) = hihc(nfl,nll) + &
        sigmah(i,nfl,nll) / brs(i,nfl,nll) / sin(ang) * &
        dels(i,nfl,nll) / 1.e5
        hidv(nfl,nll) = hidv(nfl,nll) + &
        fdpv * dels(i,nfl,nll)

!        if ( nfl == nf ) then
!          hipcp(nfl,nll)   = hipcp(nfl-1,nll)
!          hipcphi(nfl,nll) = hipcphi(nfl-1,nll)
!        endif

    enddo


! calculate collisional ion velocity
! not used; simply a diagnostic

!      do i = 1,nz
!        do ni = nion1,nion2
!          vsic(i,nfl,nll,ni) = vsid(i,nfl,nll,ni) / sumnuj(i,ni)
!        enddo
!      enddo

    return
    end subroutine update



