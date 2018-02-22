!******************************************
!******************************************

!             open_u

!******************************************
!******************************************

    subroutine open_u

! open output files (unformatted, except time.dat)

    open ( unit=70, file='time.dat'      ,form='formatted'   )
    open ( unit=71, file='deniu.dat'     ,form='unformatted' )
    open ( unit=72, file='tiu.dat'       ,form='unformatted' )
    open ( unit=73, file='vsiu.dat'      ,form='unformatted' )
    open ( unit=75, file='teu.dat'       ,form='unformatted' )
    open ( unit=78, file='vnu.dat'       ,form='unformatted' )
    open ( unit=92, file='dennu.dat'     ,form='unformatted' )
    open ( unit=93, file='hipcu.dat'     ,form='unformatted' )
    open ( unit=94, file='hihcu.dat'     ,form='unformatted' )
    open ( unit=95, file='sigmapu.dat'   ,form='unformatted' )
    open ( unit=96, file='sigmahu.dat'   ,form='unformatted' )
    open ( unit=97, file='sigmapicu.dat'   ,form='unformatted' )
    open ( unit=98, file='sigmahicu.dat'   ,form='unformatted' )
    open ( unit=711, file='deniu1.dat'     ,form='unformatted' )
    open ( unit=712, file='deniu2.dat'     ,form='unformatted' )
    open ( unit=713, file='deniu3.dat'     ,form='unformatted' )
    open ( unit=714, file='deniu4.dat'     ,form='unformatted' )
    open ( unit=715, file='deniu5.dat'     ,form='unformatted' )
    open ( unit=716, file='deniu6.dat'     ,form='unformatted' )
    open ( unit=717, file='deniu7.dat'     ,form='unformatted' )
    open ( unit=1718, file='deneu.dat'      ,form='unformatted' )
    open ( unit=811, file='tiu1.dat'     ,form='unformatted' )
    open ( unit=812, file='tiu2.dat'     ,form='unformatted' )
    open ( unit=813, file='tiu3.dat'     ,form='unformatted' )
    open ( unit=814, file='tiu4.dat'     ,form='unformatted' )
    open ( unit=815, file='tiu5.dat'     ,form='unformatted' )
    open ( unit=816, file='tiu6.dat'     ,form='unformatted' )
    open ( unit=817, file='tiu7.dat'     ,form='unformatted' )
    open ( unit=911, file='vsiu1.dat'     ,form='unformatted' )
    open ( unit=912, file='vsiu2.dat'     ,form='unformatted' )
    open ( unit=913, file='vsiu3.dat'     ,form='unformatted' )
    open ( unit=914, file='vsiu4.dat'     ,form='unformatted' )
    open ( unit=915, file='vsiu5.dat'     ,form='unformatted' )
    open ( unit=916, file='vsiu6.dat'     ,form='unformatted' )
    open ( unit=917, file='vsiu7.dat'     ,form='unformatted' )
    open ( unit=1711, file='dennu1.dat'     ,form='unformatted' )
    open ( unit=1712, file='dennu2.dat'     ,form='unformatted' )
    open ( unit=1713, file='dennu3.dat'     ,form='unformatted' )
    open ( unit=1714, file='dennu4.dat'     ,form='unformatted' )
    open ( unit=1715, file='dennu5.dat'     ,form='unformatted' )
    open ( unit=569,  file='rhsegv.dat'     ,form='unformatted' )


    open ( unit=196, file='vnqu.dat'     ,form='unformatted' )
    open ( unit=197, file='vnpu.dat'     ,form='unformatted' )
    open ( unit=198, file='vnphiu.dat'   ,form='unformatted' )
    open ( unit=201, file='jpu.dat'     ,form='unformatted' )
    open ( unit=202, file='jphiu.dat'   ,form='unformatted' )

    open ( unit=491, file='hipcpu.dat'    ,form='unformatted' )
    open ( unit=492, file='hipcphiu.dat'  ,form='unformatted' )
    open ( unit=493, file='hihcmu.dat'    ,form='unformatted' )
    open ( unit=494, file='hidpvu.dat'     ,form='unformatted' )
    open ( unit=495, file='hidphivu.dat'   ,form='unformatted' )
    open ( unit=496, file='hidpgu.dat'     ,form='unformatted' )
    open ( unit=497, file='hidphigu.dat'   ,form='unformatted' )
    open ( unit=498, file='phiu.dat'       ,form='unformatted' )

! diagnostic files (unformatted)

    open ( unit=81, file='t1u.dat'  ,form='unformatted' )
    open ( unit=82, file='t2u.dat'  ,form='unformatted' )
    open ( unit=83, file='t3u.dat'  ,form='unformatted' )
    open ( unit=84, file='u1u.dat'  ,form='unformatted' )
    open ( unit=85, file='u2u.dat'  ,form='unformatted' )
    open ( unit=86, file='u3u.dat'  ,form='unformatted' )
    open ( unit=87, file='u4u.dat'  ,form='unformatted' )
    open ( unit=88, file='u5u.dat'  ,form='unformatted' )
    open ( unit=384, file='u1pu.dat'  ,form='unformatted' )
    open ( unit=385, file='u2su.dat'  ,form='unformatted' )
    open ( unit=386, file='u3hu.dat'  ,form='unformatted' )

    return
    end subroutine open_u

!******************************************
!******************************************

!             output

!******************************************
!******************************************

    subroutine output ( hr,ntm,istep,phi,denit,dennt,vsit,sumvsit, &
    tit,ut,vt,vpit,tet,tnt,u1t, &
    u2t,u3t,u4t,vnqt,vnpt,vnphit,jpt,jphit, &
    u1pt,u2st,u3ht,sigmapict,sigmahict, &
    sigmapt,sigmaht )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use conductance_mod

    real :: denit(nz,nf,nlt,nion)
    real :: dennt(nz,nf,nlt,nion)
    real :: vsit(nz,nf,nlt,nion)
    real :: sumvsit(nz,nf,nlt,nion)
    real :: tet(nz,nf,nlt),tit(nz,nf,nlt,nion),tnt(nz,nf,nlt)
    real :: ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt)
    real :: u1t(nz,nf,nlt),u2t(nz,nf,nlt),u3t(nz,nf,nlt),u4t(nz,nf,nlt)
    real :: vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt)
    real :: jpt(nz,nf,nlt),jphit(nz,nf,nlt)
    real :: phi(nnx,nny)
    real :: u1pt(nz,nf,nlt),u2st(nz,nf,nlt),u3ht(nz,nf,nlt)
    real :: sigmapict(nz,nf,nlt),sigmahict(nz,nf,nlt)
    real :: sigmapt(nz,nf,nlt),sigmaht(nz,nf,nlt)

    real, dimension(:,:,:), allocatable :: denit1,denit2,denit3, &
    denit4,denit5,denit6,denit7,denet
    real, dimension(:,:,:), allocatable :: dennt1,dennt2,dennt3, &
    dennt4,dennt5,dennt6,dennt7
    real, dimension(:,:,:), allocatable :: tit1,tit2,tit3, &
    tit4,tit5,tit6,tit7
    real, dimension(:,:,:), allocatable ::vsit1,vsit2,vsit3, &
    vsit4,vsit5,vsit6,vsit7


    hr24   = mod (hr,24.)
    totsec = hr24 * 3600.
    thr    = totsec / 3600.
    nthr   = int(thr)
    tmin   = ( thr - nthr ) * 60.
    ntmin  = int(mod(tmin,60.))
    tsec   = ( tmin - ntmin ) * 60.
    ntsec  = int(tsec)

    print *,'istep = ',istep,' ntm = ',ntm
    print *,' hr = ',hr,' dt = ',dt

    write (70,100) ntm,nthr,ntmin,ntsec,hr

    allocate &
    (denit1(nz,nf,nlt),denit2(nz,nf,nlt),denit3(nz,nf,nlt), &
    denit4(nz,nf,nlt),denit5(nz,nf,nlt),denit6(nz,nf,nlt), &
    denit7(nz,nf,nlt),denet(nz,nf,nlt))


    do k = 1,nlt
        do j = 1,nf
            do i = 1,nz
                denit1(i,j,k) = denit(i,j,k,1)
                denit2(i,j,k) = denit(i,j,k,2)
                denit3(i,j,k) = denit(i,j,k,3)
                denit4(i,j,k) = denit(i,j,k,4)
                denit5(i,j,k) = denit(i,j,k,5)
                denit6(i,j,k) = denit(i,j,k,6)
                denit7(i,j,k) = denit(i,j,k,7)
            enddo
        enddo
    enddo

    open (144,file='denit_cg.rst',form='unformatted')
    write(144) denit1,denit2,denit3,denit4,denit5,denit6,denit7
    close(144)

    do k = 1,nlt
        do j = 1,nf
            do i = 1,nz
                denet(i,j,k) = 0
                do ni = nion1,nion2
                    denet(i,j,k) = denet(i,j,k) + denit(i,j,k,ni)
                enddo
            enddo
        enddo
    enddo

        write(1718) denet
        write(711) denit1
        write(712) denit2
        write(713) denit3
        write(714) denit4
        write(715) denit5
    !         write(716) denit6
    !         write(717) denit7

    deallocate (denit1,denit2,denit3, &
                denit4,denit5,denit6,denit7,denet)

    allocate &
    (dennt1(nz,nf,nlt),dennt2(nz,nf,nlt),dennt3(nz,nf,nlt), &
    dennt4(nz,nf,nlt),dennt5(nz,nf,nlt),dennt6(nz,nf,nlt), &
    dennt7(nz,nf,nlt))

    do k = 1,nlt
        do j = 1,nf
            do i = 1,nz
                dennt1(i,j,k) = dennt(i,j,k,1)
                dennt2(i,j,k) = dennt(i,j,k,2)
                dennt3(i,j,k) = dennt(i,j,k,3)
                dennt4(i,j,k) = dennt(i,j,k,4)
                dennt5(i,j,k) = dennt(i,j,k,5)
                dennt6(i,j,k) = dennt(i,j,k,6)
                dennt7(i,j,k) = dennt(i,j,k,7)
            enddo
        enddo
    enddo


        write(1711) dennt1
        write(1712) dennt2
        write(1713) dennt3
        write(1714) dennt4
        write(1715) dennt5

    deallocate (dennt1,dennt2,dennt3, &
    dennt4,dennt5,dennt6,dennt7)

    allocate &
    (tit1(nz,nf,nlt),tit2(nz,nf,nlt),tit3(nz,nf,nlt), &
    tit4(nz,nf,nlt),tit5(nz,nf,nlt),tit6(nz,nf,nlt), &
    tit7(nz,nf,nlt))

    do k = 1,nlt
        do j = 1,nf
            do i = 1,nz
                tit1(i,j,k) = tit(i,j,k,1)
                tit2(i,j,k) = tit(i,j,k,2)
                tit3(i,j,k) = tit(i,j,k,3)
                tit4(i,j,k) = tit(i,j,k,4)
                tit5(i,j,k) = tit(i,j,k,5)
                tit6(i,j,k) = tit(i,j,k,6)
                tit7(i,j,k) = tit(i,j,k,7)
            enddo
        enddo
    enddo

        write(811) tit1
        write(812) tit2
        write(813) tit3
    !         write(814) tit4
    !         write(815) tit5
    !         write(816) tit6
    !         write(817) tit7

        write(75) tet

    deallocate (tit1,tit2,tit3,tit4,tit5,tit6,tit7)

    allocate &
    (vsit1(nz,nf,nlt),vsit2(nz,nf,nlt),vsit3(nz,nf,nlt), &
    vsit4(nz,nf,nlt),vsit5(nz,nf,nlt),vsit6(nz,nf,nlt), &
    vsit7(nz,nf,nlt))


    do k = 1,nlt
        do j = 1,nf
            do i = 1,nz
                vsit1(i,j,k) = vsit(i,j,k,1)
                vsit2(i,j,k) = vsit(i,j,k,2)
                vsit3(i,j,k) = vsit(i,j,k,3)
                vsit4(i,j,k) = vsit(i,j,k,4)
                vsit5(i,j,k) = vsit(i,j,k,5)
                vsit6(i,j,k) = vsit(i,j,k,6)
                vsit7(i,j,k) = vsit(i,j,k,7)
            enddo
        enddo
    enddo

        write(911) vsit1
        write(912) vsit2
        write(913) vsit3
    !         write(914) vsit4
        write(915) vsit5
    !         write(916) vsit6
    !         write(917) vsit7


    deallocate (vsit1,vsit2,vsit3, &
    vsit4,vsit5,vsit6,vsit7)

    !         write(78) vnt
    !         write(81) t1t
    !         write(82) t2t
    !         write(83) t3t
        write(384) u1pt
        write(385) u2st
        write(386) u3ht
        write(84) u1t
        write(85) u2t
        write(86) u3t
        write(87) u4t
    !         write(88) u5t
    !         write(93) hipct
    !         write(94) hihct
        write(95) sigmapt
        write(96) sigmaht
        write(97) sigmapict
        write(98) sigmahict

        write(196) vnqt
        write(197) vnpt
        write(198) vnphit
        write(201) jpt
        write(202) jphit

        write(491) hipcpt
        write(492) hipcphit
        write(493) hihcmt
        write(494) hidpvt
        write(495) hidphivt
        write(496) hidpgt
        write(497) hidphigt
        write(498) phi

    !         write(569) cxxe,cyye,cxe,cye,rhse
    !         write(569) rhseg,rhsev
    !        close(569)

    100 format(1x,4i6,1p1e14.4)
    101 format(1x,1p10e16.6)

    return
    end subroutine output

