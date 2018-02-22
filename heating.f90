!******************************************
!******************************************
     
!            htemp
     
!******************************************
!******************************************
     
    subroutine htemp ( tti,tiold,tvn,nuin,nfl,nll )
     
    use parameter_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use exb_mod
    use grid_mod

    real :: tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
    real :: tvn(nz,nl),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
    real :: lambda
    real :: divvexb(nz)
     
    convfac = amu / bolt / 3.
     
    do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
    enddo

    do i = 1,nz
         
    ! from schunk/nagy book
         
        lambda = 23. - &
        &            0.5*alog(deni(i,nfl,nll,pth)/ &
        (ti(i,nfl,nll,pth)/evtok)**3)
        schunkfac = 0.
        do ni = nion1,nion2
            if (ni /= pth) then
                schunkfac = schunkfac + &
                deni(i,nfl,nll,ni)/deni(i,nfl,nll,pth) * &
                sqrt(ami(ni)/(ami(ni)+ami(pth))**5) * &
                (3*ami(pth)**2 + 1.6*ami(pth)*ami(ni) + 1.3*ami(ni)**2)
            endif
        enddo
        kapi(i) = 15.*3.1e4 * sqrt ( ti(i,nfl,nll,pth)**5 ) / &
        sqrt(ami(pth)) / &
        (1. + 1.75*schunkfac) / lambda
         
        kapi(i)  = 0.6667 * kapi(i) * evtok
         
    ! neutrals
         
        do nn = 1,nneut
            redmass = &
            ami(pthp) * amn(nn) / ( ami(pthp) + amn(nn) ) ** 2
            s2i(i) = s2i(i) + 2. * nuin(i,pthp,nn) * redmass
            s3i(i) = s3i(i) &
            + convfac * amn(nn) &
            * abs ( vsi(i,nfl,nll,pthp) - tvn(i,nll) ) ** 2 &
            * 2. * nuin(i,pthp,nn) * redmass
        enddo
         
        s1i(i) = s2i(i) * tn(i,nfl,nll)
         
    ! electrons
         
        s4i(i) = 7.7e-6 * ne(i,nfl,nll) / ami(pthp) &
        / te(i,nfl,nll) / sqrt(te(i,nfl,nll)) &
        * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl,nll)
         
    ! other ions
         
        do ni = nion1,nion2
            if ( ni /= pthp ) then
                tfac    =    ti(i,nfl,nll,pthp) / ami(pthp) &
                +  ti(i,nfl,nll,ni) / ami(ni)
                xs6i    = 3.3e-4 * deni(i,nfl,nll,ni) / ami(pthp) / ami(ni) &
                / tfac / sqrt(tfac) * .66667 * evtok
                xs7i    = xs6i * ti(i,nfl,nll,ni)
                s6i(i) = s6i(i) + xs6i
                s7i(i) = s7i(i) + xs7i
            endif
        enddo
         
    enddo
     
! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift
     
    nzh = nz / 2
    vexbeq = vexbp(nzh,nfl,nll)
    do i = 1,nz
        divvexb(i) = 6.*vexbeq / &
        (ps(i,nfl,nll)*re*1.e5) * &
        cos(blats(i,nfl,nll)*po180)**2 * &
        (1.+sin(blats(i,nfl,nll)*po180)**2) / &
        (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2
        s2i(i) = s2i(i) - 0.3333 * divvexb(i)
    enddo
     
    call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthp, &
    nfl,nll)
     
    return
    end subroutine htemp


!******************************************
!******************************************

!            hetemp

!******************************************
!******************************************

    subroutine hetemp ( tti,tiold,tvn,nuin,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use exb_mod
    use grid_mod


    real :: tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
    real :: tvn(nz,nl),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
    real :: lambda
    real :: divvexb(nz)

    convfac = amu / bolt / 3.

    do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
    enddo

    do i = 1,nz

    ! from schunk/nagy book

        lambda = 23. - &
        &            0.5*alog(deni(i,nfl,nll,pthe)/ &
        (ti(i,nfl,nll,pthe)/evtok)**3)
        schunkfac = 0.
        do ni = nion1,nion2
            if (ni /= pthe) then
                schunkfac = schunkfac + &
                deni(i,nfl,nll,ni)/deni(i,nfl,nll,pthe) * &
                sqrt(ami(ni)/(ami(ni)+ami(pthe))**5) * &
                (3*ami(pthe)**2 + 1.6*ami(pthe)*ami(ni) + &
                &              1.3*ami(ni)**2)
            endif
        enddo
        kapi(i) = 15.*3.1e4 * sqrt ( ti(i,nfl,nll,pthe)**5 ) / &
        sqrt(ami(pthe)) / &
        (1. + 1.75*schunkfac) / lambda

        kapi(i)  = 0.6667 * kapi(i) * evtok

    ! neutrals

        do nn = 1,nneut
            redmass = &
            ami(pthep) * amn(nn) / ( ami(pthep) + amn(nn) ) ** 2
            s2i(i) = s2i(i) + 2. * nuin(i,pthep,nn) * redmass
            s3i(i) = s3i(i) &
            + convfac * amn(nn) &
            * abs ( vsi(i,nfl,nll,pthep) - tvn(i,nll) ) ** 2 &
            * 2. * nuin(i,pthep,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl,nll)

    ! electrons

        s4i(i) = 7.7e-6 * ne(i,nfl,nll) / ami(pthep) &
        / te(i,nfl,nll) / sqrt(te(i,nfl,nll)) &
        * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl,nll)

    ! other ions

        do ni = nion1,nion2
            if ( ni /= pthep ) then
                tfac    =   ti(i,nfl,nll,pthep) / ami(pthep) &
                + ti(i,nfl,nll,ni) / ami(ni)
                xs6i    = 3.3e-4 * deni(i,nfl,nll,ni) / ami(pthep) / ami(ni) &
                / tfac / sqrt(tfac) * .66667 * evtok
                xs7i    = xs6i * ti(i,nfl,nll,ni)
                s6i(i) = s6i(i) + xs6i
                s7i(i) = s7i(i) + xs7i
            endif
        enddo

    enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

    nzh = nz / 2
    vexbeq = vexbp(nzh,nfl,nll)
    do i = 1,nz
        divvexb(i) = 6.*vexbeq / &
        (ps(i,nfl,nll)*re*1.e5) * &
        cos(blats(i,nfl,nll)*po180)**2 * &
        (1.+sin(blats(i,nfl,nll)*po180)**2) / &
        (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2
        s2i(i) = s2i(i) - 0.3333 * divvexb(i)
    enddo

    call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthep, &
    nfl,nll)

    return
    end subroutine hetemp

!******************************************
!******************************************

!            otemp

!******************************************
!******************************************

    subroutine otemp ( tti,tiold,tvn,nuin,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use exb_mod
    use grid_mod

    real :: tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
    real :: tvn(nz,nl),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
    real :: lambda
    real :: divvexb(nz)

    convfac = amu / bolt / 3.

    do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
    enddo

    do i = 1,nz

    ! from schunk/nagy book

        lambda = 23. - &
        &            0.5*alog(deni(i,nfl,nll,pto)/ &
        (ti(i,nfl,nll,pto)/evtok)**3)
        schunkfac = 0.
        do ni = nion1,nion2
            if (ni /= pto) then
                schunkfac = schunkfac + &
                deni(i,nfl,nll,ni)/deni(i,nfl,nll,pto) * &
                sqrt(ami(ni)/(ami(ni)+ami(pto))**5) * &
                (3*ami(pto)**2 + 1.6*ami(pto)*ami(ni) + 1.3*ami(ni)**2)
            endif
        enddo
        kapi(i) = 15.*3.1e4 * sqrt ( ti(i,nfl,nll,pto)**5 ) / &
        sqrt(ami(pto)) / &
        (1. + 1.75*schunkfac) / lambda

        kapi(i)  = 0.6667 * kapi(i) * evtok

    ! neutrals

        do nn = 1,nneut
            redmass = &
            ami(ptop) * amn(nn) / ( ami(ptop) + amn(nn) ) ** 2
            s2i(i) = s2i(i) + 2. * nuin(i,ptop,nn) * redmass
            s3i(i) = s3i(i) &
            + convfac * amn(nn) &
            * abs ( vsi(i,nfl,nll,ptop) - tvn(i,nll) ) ** 2 &
            * 2. * nuin(i,ptop,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl,nll)

    ! electrons

        s4i(i) = 7.7e-6 * ne(i,nfl,nll) / ami(ptop) &
        / te(i,nfl,nll) / sqrt(te(i,nfl,nll)) &
        * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl,nll)

    ! other ions

        do ni = nion1,nion2
            if ( ni /= ptop ) then
                tfac    =    ti(i,nfl,nll,ptop) / ami(ptop) &
                + ti(i,nfl,nll,ni) / ami(ni)
                xs6i    = 3.3e-4 * deni(i,nfl,nll,ni) / ami(ptop) / ami(ni) &
                / tfac / sqrt(tfac) * .66667 * evtok
                xs7i    = xs6i * ti(i,nfl,nll,ni)
                s6i(i) = s6i(i) + xs6i
                s7i(i) = s7i(i) + xs7i
            endif
        enddo

    enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

    nzh = nz / 2
    vexbeq = vexbp(nzh,nfl,nll)
    do i = 1,nz
        divvexb(i) = 6.*vexbeq / &
        (ps(i,nfl,nll)*re*1.e5) * &
        cos(blats(i,nfl,nll)*po180)**2 * &
        (1.+sin(blats(i,nfl,nll)*po180)**2) / &
        (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2
        s2i(i) = s2i(i) - 0.3333 * divvexb(i)
    enddo

    call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,ptop, &
    nfl,nll)

    return
    end subroutine otemp



!******************************************
!******************************************

!            etemp

!******************************************
!******************************************

    subroutine etemp ( tte,te_old,phprodr,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use exb_mod
    use grid_mod

    real :: tte(nz),te_old(nz),kape(nz)
    real :: s1e(nz),s2e(nz),s3e(nz),s4e(nz),phprodr(nz,nion)
    real :: s5e(nz),qphe(nz),phprod(nz)
    real :: qen(nz,nneut)
    real :: ne300s,ne300n,n2300
    real :: ratio(nz)
    real :: divvexb(nz)
    integer :: iz300s(nf,nl),iz300n(nf,nl)

    do i = 1,nz
        s1e(i)  = 0.
        s2e(i)  = 0.
        s3e(i)  = 0.
        s4e(i)  = 0.
        kape(i) = 0.
        do ni = 1,nneut
            qen(i,ni) = 0.
        enddo
    enddo

    do i = 1,nz

        fac1 = denn(i,nfl,nll,pto)  * 1.1e-16 &
        * ( 1. + 5.7e-4 * te(i,nfl,nll) )
        fac2 = denn(i,nfl,nll,ptn2) * 2.82e-17 &
        * ( 1  - 1.2e-4 * te(i,nfl,nll) )* sqrt(te(i,nfl,nll))
        fac3 = denn(i,nfl,nll,pto2) * 2.2e-16 &
        * ( 1. + 3.6e-2  * sqrt(te(i,nfl,nll)) )
        akpefac = fac1 + fac2 + fac3

        kape(i) = 7.7e5 * sqrt ( te(i,nfl,nll)**5 ) * 0.6667 * evtok &
        / ( 1. + 3.22e4 * ( te(i,nfl,nll)**2 / &
        ne(i,nfl,nll) * akpefac) )


    ! neutrals (Tn - Te) term

    ! N2

    ! vibrational state from red book (p. 269) milward et al.

        qen(i,ptn2) = .6667 *  evtok * denn(i,nfl,nll,ptn2) * &
        ( 1.2e-19 * ( 1. - 1.2e-4 * te(i,nfl,nll) ) &
        * te(i,nfl,nll) + &
        &                     2.e-14 / sqrt(te(i,nfl,nll)) &
        + 6.5e-22 * ( tn(i,nfl,nll) - 310 ) ** 2 * &
        exp(.0023*(te(i,nfl,nll) - tn(i,nfl,nll))))

    ! O2

        qen(i,pto2) = .6667 * evtok * denn(i,nfl,nll,pto2) * &
        ( 7.9e-19 * ( 1. + 3.6e-2 * sqrt(te(i,nfl,nll))) &
        *  sqrt(te(i,nfl,nll)) + &
        &                    7.e-14 / sqrt(te(i,nfl,nll)) )

    ! O

        qen(i,pto) = .6667 * 7.2e-18 * evtok * denn(i,nfl,nll,pto) * &
        sqrt(te(i,nfl,nll))

    ! H

        qen(i,pth) = .6667 * 6.3e-16 * evtok * denn(i,nfl,nll,pth) * &
        ( 1. - 1.35e-4 * te(i,nfl,nll) ) * &
        sqrt(te(i,nfl,nll))

        do nn = 1,nneut
            s2e(i) = s2e(i) + qen(i,nn)
        enddo

        s1e(i) = s2e(i) * tn(i,nfl,nll)

    ! ions (Ti - Te) term

        do ni = nion1,nion2
            xs3e    = 7.7e-6 * deni(i,nfl,nll,ni) / ami(ni) &
            / te(i,nfl,nll) / sqrt(te(i,nfl,nll)) &
            * .66667 * evtok
            xs4e    = xs3e * ti(i,nfl,nll,ni)
            s3e(i) = s3e(i) + xs3e
            s4e(i) = s4e(i) + xs4e
        enddo

    enddo

! photoelectron heating
! red book (millward et al. p. 269)

! calculate total ion photoproduction (= photoelectron)

    do i = 1,nz
        phprod(i)   = 0.
        do ni = nion1,nion2
            phprod(i) = phprodr(i,ni) * denn(i,nfl,nll,ni) + phprod(i)
        enddo
    enddo

! iz300s/iz300n are redefined here

    do i = 1,nz
        ratio(i) = ne(i,nfl,nll) / &
        (0.1*denn(i,nfl,nll,pto)+ &
        denn(i,nfl,nll,pto2)+denn(i,nfl,nll,ptn2))
    enddo

    i = 1
    do while ( ratio(i) <= 3.e-3 .AND. i < nz )
        iz300s(nfl,nll) = i
        i         = i + 1
    enddo
     
    i = nz
    do while ( ratio(i) <= 3.e-3 .AND. i > 1 )
        iz300n(nfl,nll) = i
        i         = i - 1
    enddo

    if ( iz300s(nfl,nll) > iz300n(nfl,nll) ) then

        do i = 1,nz
            xarg =   ne(i,nfl,nll) &
            / (        denn(i,nfl,nll,pto2) &
            +      denn(i,nfl,nll,ptn2) &
            + .1 * denn(i,nfl,nll,pto)   )
            x    = alog ( xarg )
            earg =     12.75 &
            + 6.941 * x &
            + 1.166 * x ** 2 &
            + 0.08034 * x ** 3 &
            + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
        enddo
    else
        do i = 1,iz300s(nfl,nll)
            xarg =   ne(i,nfl,nll) &
            / (        denn(i,nfl,nll,pto2) &
            +      denn(i,nfl,nll,ptn2) &
            + .1 * denn(i,nfl,nll,pto)   )
            x    = alog ( xarg )
            earg =     12.75 &
            + 6.941 * x &
            + 1.166 * x ** 2 &
            + 0.08034 * x ** 3 &
            + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
        enddo

    ! smooth things at '300' km

        izs   = iz300s(nfl,nll)
        facts = (3.e-3-ratio(izs)) / &
        (ratio(izs+1)-ratio(izs))
        ne300s = ne(izs,nfl,nll) + (ne(izs+1,nfl,nll)- &
        ne(izs,nfl,nll)) * facts
        o2300 = denn(izs,nfl,nll,pto2) + &
        (denn(izs+1,nfl,nll,pto2)- &
        denn(izs,nfl,nll,pto2)) * facts
        n2300 = denn(izs,nfl,nll,ptn2) + &
        (denn(izs+1,nfl,nll,ptn2)- &
        denn(izs,nfl,nll,ptn2)) * facts
        o300 = denn(izs,nfl,nll,pto) + &
        (denn(izs+1,nfl,nll,pto)-denn(izs,nfl,nll,pto)) * facts
        phprod300 = phprod(izs) + &
        (phprod(izs+1)-phprod(izs)) * facts
        xarg300 = ne300s / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 + &
        &         6.941 * x300 + &
        &         1.166 * x300 ** 2 + &
        &         0.08034 * x300 ** 3 + &
        &         0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0s = epsi300 * phprod300 / ne300s

        do i = iz300n(nfl,nll),nz
            xarg =   ne(i,nfl,nll) &
            / (       denn(i,nfl,nll,pto2) &
            +      denn(i,nfl,nll,ptn2) &
            + .1 * denn(i,nfl,nll,pto) )
            x    = alog ( xarg )
            earg =     12.75 &
            + 6.941 * x &
            + 1.166 * x ** 2 &
            + 0.08034 * x ** 3 &
            + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
        enddo

        izn   = iz300n(nfl,nll)
        factn = (3.e-3-ratio(izn)) / &
        (ratio(izn-1)-ratio(izn))
        ne300n = ne(izn,nfl,nll) + &
        (ne(izn-1,nfl,nll)-ne(izn,nfl,nll)) * factn
        o2300 = denn(izn,nfl,nll,pto2) + &
        (denn(izn-1,nfl,nll,pto2)- &
        denn(izn,nfl,nll,pto2)) * factn
        n2300 = denn(izn,nfl,nll,ptn2) + &
        (denn(izn-1,nfl,nll,ptn2)- &
        denn(izn,nfl,nll,ptn2)) * factn
        o300 = denn(izn,nfl,nll,pto) + &
        (denn(izn-1,nfl,nll,pto)-denn(izn,nfl,nll,pto)) * factn
        phprod300 = phprod(izn) + &
        (phprod(izn-1)-phprod(izn)) * factn
        xarg300 = ne300n / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 + &
        &         6.941 * x300 + &
        &         1.166 * x300 ** 2 + &
        &         0.08034 * x300 ** 3 + &
        &         0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0n = epsi300 * phprod300 / ne300n

        xbms = bms(izs,nfl,nll) + &
        (bms(izs+1,nfl,nll)-bms(izs,nfl,nll)) * facts
        xbmn = bms(izn,nfl,nll) + &
        (bms(izn-1,nfl,nll)-bms(izn,nfl,nll)) * factn

        dels300s = dels(iz300s(nfl,nll),nfl,nll) * facts
        dels300n = dels(iz300n(nfl,nll)-1,nfl,nll) * factn

    ! MS: Old code used a wasteful way to calculate xn.
    ! Cleaner version here.
        xn = 0.
    ! Set bottom integration bound to 300 km.
        xn =   xn + 0.5 * ( ne(iz300n(nfl,nll)-1,nfl,nll) + ne300n ) * &
        (dels(iz300n(nfl,nll)-1,nfl,nll) - dels300n )
        do i =iz300n(nfl,nll)-2,iz300s(nfl,nll)+1,-1
            xn = xn + 0.5 * ( ne(i,nfl,nll) + ne(i+1,nfl,nll) ) * &
            dels(i,nfl,nll)
        enddo

        if ( q0s < 0 .OR. q0n < 0 ) then
            print *,' q0s = ',q0s,' q0n = ',q0n,' nfl = ',nfl
        endif

    ! 1/22/00

    ! put in dels (arc length along field line)

        xs    = 0.
        do i = iz300s(nfl,nll)+1,iz300n(nfl,nll)-1
            if (i == iz300s(nfl,nll)+1) then
                xs = xs + 0.5*( ne300s + ne(i,nfl,nll) ) * &
                (dels(iz300s(nfl,nll),nfl,nll) - dels300s)
            else
                xs = xs + 0.5 * ( ne(i,nfl,nll) + ne(i-1,nfl,nll) ) &
                * dels(i-1,nfl,nll)
                xn = xn - 0.5 * ( ne(i,nfl,nll) + ne(i-1,nfl,nll) ) &
                * dels(i-1,nfl,nll)
            endif
             
            xints = cqe*xs
            xintn = cqe*xn
            xqs    = ne(i,nfl,nll) * q0s * bms(i,nfl,nll) / &
            xbms * exp(-xints)
            xqn    = ne(i,nfl,nll) * q0n * bms(i,nfl,nll) / &
            xbmn * exp(-xintn)
            qphe(i) = xqs + xqn
        enddo
    endif

    do i = 1,nz
        s5e(i) = 0.66667 * evtok * qphe(i) / ne(i,nfl,nll) !* .15
    enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

    nzh    = nz / 2
    vexbeq = vexbp(nzh,nfl,nll)
    do i = 1,nz
        divvexb(i) = 6.*vexbeq / &
        (ps(i,nfl,nll)*re*1.e5) * &
        cos(blats(i,nfl,nll)*po180)**2 * &
        (1.+sin(blats(i,nfl,nll)*po180)**2) / &
        (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2
        s2e(i) = s2e(i) - 0.3333 * divvexb(i)
    enddo

    call tesolv(tte,te_old,kape,s1e,s2e,s3e,s4e,s5e,nfl,nll)

    return
    end subroutine etemp

