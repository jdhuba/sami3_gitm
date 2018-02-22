!******************************************
!******************************************

!            neutambt

!******************************************
!******************************************


! calculate neutral densities and temperature
! from nrlmsise00

! no obtained from eq. (128) - bailey and balan (red book)

! neutral density and temperature

! input:
!    iyd - year and day as yyddd
!    sec - ut(sec)
!    alt - altitude(km) (greater than 85 km)
!    glat - geodetic latitude(deg)
!    glong - geodetic longitude(deg)
!    stl - local apparent solar time(hrs)
!    f107a - 3 month average of f10.7 flux
!    f107 - daily f10.7 flux for previous day
!    ap - magnetic index(daily) or when sw(9)=-1. :
!       - array containing:
!         (1) daily ap
!         (2) 3 hr ap index for current time
!         (3) 3 hr ap index for 3 hrs before current time
!         (4) 3 hr ap index for 6 hrs before current time
!         (5) 3 hr ap index for 9 hrs before current time
!         (6) average of eight 3 hr ap indicies from 12 to 33 hrs prior
!             to current time
!         (7) average of eight 3 hr ap indicies from 36 to 59 hrs prior
!             to current time
!    mass - mass number (only density for selected gas is
!             calculated.  mass 0 is temperature.  mass 48 for all.
! output:
!    d(1) - he number density(cm-3)
!    d(2) - o number density(cm-3)
!    d(3) - n2 number density(cm-3)
!    d(4) - o2 number density(cm-3)
!    d(5) - ar number density(cm-3)
!    d(6) - total mass density(gm/cm3)
!    d(7) - h number density(cm-3)
!    d(8) - n number density(cm-3)
!    d(9) - anomalous O (see msis)
!    t(1) - exospheric temperature
!    t(2) - temperature at alt

! neutral winds

!    iyd - year and day as yyddd
!    sec - ut(sec)  (not important in lower atmosphere)
!    alt - altitude(km)
!    glat - geodetic latitude(deg)
!    glong - geodetic longitude(deg)
!    stl - local apparent solar time(hrs)
!    f107a - 3 month average of f10.7 flux (use 150 in lower atmos.)
!    f107 - daily f10.7 flux for previous day ( " )
!    ap - two element array with
!         ap(1) = magnetic index(daily) (use 4 in lower atmos.)
!         ap(2)=current 3hr ap index (used only when sw(9)=-1.)
! note:  ut, local time, and longitude are used independently in the
!        model and are not of equal importance for every situation.
!        for the most physically realistic calculation these three
!        variables should be consistent.
! output
!    w(1) = meridional (m/sec + northward)
!    w(2) = zonal (m/sec + eastward)


    subroutine neutambt (hr,nll)

    use gitm_mod
    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use atomic_mod
    use grid_mod
    use message_passing_mod

    real :: d(9),temp(2)
    real :: whm93(2),app(2)

    hruti = hr

! GITM

      hruti24 = mod(hruti,24.)

    do j = 1,nf
        do i = 1,nz
            glonsij = glons(i,j,nll)
            if ( lcr ) glonsij = glons0(i,j,nll) - (hruti - hrinit) * 15.
            glonsij = mod(glonsij,360.)
            hrl   = mod(hruti + glonsij / 15.,24.)
            call msistim ( int(year),int(day),hrl, &
                           glonsij,iyd,sec )
            call gtd7 ( iyd,sec,alts(i,j,nll), &
                        glats(i,j,nll),glonsij, &
                        hrl,fbar,f10p7,aap,mmass,d,temp )
            denni(i,j,nll,pth )  = snn(pth)  * d(7)
            denni(i,j,nll,pthe)  = snn(pthe) * d(1)
            denni(i,j,nll,ptn )  = snn(ptn)  * d(8)
            denni(i,j,nll,pto )  = snn(pto)  * d(2)
            denni(i,j,nll,ptn2)  = snn(ptn2) * d(3) + 1.e-30
            denni(i,j,nll,pto2)  = snn(pto2) * d(4) + 1.e-30
            tni(i,j,nll)         = stn * temp(2)
            denni(i,j,nll,ptno)  = 0.4 * exp( -3700. / tni(i,j,nll) ) &
                                 * denni(i,j,nll,pto2) &
                                 + 5.0e-7 * denni(i,j,nll,pto)
        enddo
    enddo
    do j = 1,nf
        do i = 1,nz
            app(1)   = ap
            app(2)   = ap
            glonsij = glons(i,j,nll)
            if ( lcr ) glonsij = glons0(i,j,nll) - (hruti - hrinit) * 15.
            glonsij = mod(glonsij,360.)
            hrl     = mod(hruti + glonsij / 15.,24.)
            call msistim ( int(year),int(day),hrl,glonsij,iyd,sec )
            if(lhwm93)  call gws5 ( iyd,sec,alts(i,j,nll), &
                                    glats(i,j,nll),glonsij, &
                                    hrl,fbar,f10p7,app,whm93        )
            if(lhwm14)  call hwm14 ( iyd,sec,alts(i,j,nll), &
                                     glats(i,j,nll),glonsij, &
                                     hrl,fbar,f10p7,app,whm93        )

            vi(i,j,nll)   = 100. * whm93(1) * tvn0 ! convert to cm/sec
            ui(i,j,nll)   = 100. * whm93(2) * tvn0 ! convert to cm/sec
            wi(i,j,nll)   = vw   * tvn0

! use GITM data

          if (alts(i,j,nll) .gt. 101. ) then

             call get_gitm_data_i(hruti24,glats(i,j,nll),glonsij,&
                         alts(i,j,nll),un1,vn1,&
!                         ttn1,tn21,to21,to11,tno1,tn4s1)
                         ttn1,tn21,to21,to11,tno1,tn4s1,thyd1)

!c convert to cm/s
 
              vi(i,j,nll) = vn1*100.
              ui(i,j,nll) = un1*100.

!c overwrite
! need to comment out these lines to use MSIS
! upper limit 2000 so hydrogen transitions to MSIS in the geocorona

            if ( alts(i,j,nll) <= 2000. ) denni(i,j,nll,pth) = thyd1
            denni(i,j,nll,ptn )  = tn4s1
            denni(i,j,nll,pto )  = to11
            denni(i,j,nll,ptn2)  = tn21
            denni(i,j,nll,pto2)  = to21
            tni(i,j,nll)         = ttn1
            denni(i,j,nll,ptno)  = tno1
          endif

        enddo
    enddo



    hrutf = hr + dt_gitm

! GITM

      hrutf24 = mod(hrutf,24.)

    do j = 1,nf
        do i = 1,nz
            glonsij = glons(i,j,nll)
            if ( lcr ) glonsij = glons0(i,j,nll) - (hrutf - hrinit) * 15.
            glonsij = mod(glonsij,360.)
            hrl   = mod(hrutf + glonsij / 15.,24.)
            call msistim ( int(year),int(day),hrl, &
                           glonsij,iyd,sec )
            call gtd7 ( iyd,sec,alts(i,j,nll), &
                        glats(i,j,nll),glonsij, &
                        hrl,fbar,f10p7,aap,mmass,d,temp )
            dennf(i,j,nll,pth )  = snn(pth)  * d(7)
            dennf(i,j,nll,pthe)  = snn(pthe) * d(1)
            dennf(i,j,nll,ptn )  = snn(ptn)  * d(8)
            dennf(i,j,nll,pto )  = snn(pto)  * d(2)
            dennf(i,j,nll,ptn2)  = snn(ptn2) * d(3) + 1.e-30
            dennf(i,j,nll,pto2)  = snn(pto2) * d(4) + 1.e-30
            tnf(i,j,nll)         = stn * temp(2)
            dennf(i,j,nll,ptno)  = 0.4 * exp( -3700. / tnf(i,j,nll) ) &
                                 * dennf(i,j,nll,pto2) &
                                 + 5.0e-7 * dennf(i,j,nll,pto)
        enddo
    enddo
    do j = 1,nf
        do i = 1,nz
            app(1)   = ap
            app(2)   = ap
            glonsij = glons(i,j,nll)
            if ( lcr ) glonsij = glons0(i,j,nll) - (hrutf - hrinit) * 15.
            glonsij = mod(glonsij,360.)
            hrl     = mod(hrutf + glonsij / 15.,24.)
            call msistim ( int(year),int(day),hrl,glonsij,iyd,sec )
            if(lhwm93)  call gws5 ( iyd,sec,alts(i,j,nll), &
                                    glats(i,j,nll),glonsij, &
                                    hrl,fbar,f10p7,app,whm93        )
            if(lhwm14)  call hwm14 ( iyd,sec,alts(i,j,nll), &
                                     glats(i,j,nll),glonsij, &
                                     hrl,fbar,f10p7,app,whm93        )

            vf(i,j,nll)   = 100. * whm93(1) * tvn0 ! convert to cm/sec
            uf(i,j,nll)   = 100. * whm93(2) * tvn0 ! convert to cm/sec
            wf(i,j,nll)   = vw   * tvn0

! use GITM data

          if (alts(i,j,nll) .gt. 101. ) then
              
             call get_gitm_data_f(hrutf24,glats(i,j,nll),glonsij,&
                         alts(i,j,nll),un1,vn1,&
!                         ttn1,tn21,to21,to11,tno1,tn4s1)
                         ttn1,tn21,to21,to11,tno1,tn4s1,thyd1)

!c convert to cm/s
 
              vf(i,j,nll) = vn1*100.
              uf(i,j,nll) = un1*100.

!c overwrite
! need to comment out these lines to use MSIS
! upper limit 2000 so hydrogen transitions to MSIS in the geocorona

            if ( alts(i,j,nll) <= 2000. ) dennf(i,j,nll,pth) = thyd1
            dennf(i,j,nll,ptn )  = tn4s1
            dennf(i,j,nll,pto )  = to11
            dennf(i,j,nll,ptn2)  = tn21
            dennf(i,j,nll,pto2)  = to21
            tnf(i,j,nll)         = ttn1
            dennf(i,j,nll,ptno)  = tno1
          endif

        enddo
    enddo

!     set density, temperature, velocity to current time

    do k = 1,nneut
        do n = 1,nl
            do j = 1,nf
                do i = 1,nz
                    denn(i,j,n,k)  = denni(i,j,n,k)
                enddo
            enddo
        enddo
    enddo

    do n = 1,nl
        do j = 1,nf
            do i = 1,nz
                tn(i,j,n)      = tni(i,j,n)
                u(i,j,n)       = ui(i,j,n)
                v(i,j,n)       = vi(i,j,n)
                w(i,j,n)       = wi(i,j,n)
            enddo
        enddo
    enddo

!    print *,'some u,v',taskid,u(nz/2,10,nl/2),v(nz/2,10,nl/2)


!!$    if ( taskid == 1 ) print *,'some max ',maxval(u),maxval(v),&
!!$            maxval(denn(:,:,:,pth)),maxval(denn(:,:,:,ptn)),&
!!$            maxval(denn(:,:,:,pto)),maxval(denn(:,:,:,ptn2)),&
!!$            maxval(denn(:,:,:,pto2)),maxval(denn(:,:,:,ptno)),&
!!$            maxval(tn)

    return
    end subroutine neutambt





!******************************************
!******************************************

!            neut

!******************************************
!******************************************

    subroutine neut(hrut)

    use parameter_mod
    use variable_mod
    use time_mod

    tfactor  = ( hrut - hruti ) / ( hrutf - hruti )
    tfactor1 = 1. - tfactor

    do k = 1,nneut
        do n = 1,nl
            do j = 1,nf
                do i = 1,nz
                    denn(i,j,n,k) = denni(i,j,n,k) * tfactor1 &
                                  + dennf(i,j,n,k) * tfactor
                enddo
            enddo
        enddo
    enddo

    do n = 1,nl
        do j = 1,nf
            do i = 1,nz
                tn(i,j,n) =   tni(i,j,n) * tfactor1 &
                            + tnf(i,j,n) * tfactor
                u(i,j,n)  =   ui(i,j,n)  * tfactor1 &
                            + uf(i,j,n)  * tfactor
                v(i,j,n)  =   vi(i,j,n)  * tfactor1 &
                            + vf(i,j,n)  * tfactor
                w(i,j,n)  =   wi(i,j,n)  * tfactor1 &
                            + wf(i,j,n)  * tfactor
            enddo
        enddo
    enddo

    return
    end subroutine neut




!******************************************
!******************************************

!            msistim

!******************************************
!******************************************

    subroutine msistim ( iyr,iday,hr,glong,iyd,secut )

! msistim calculates time parameters for the
! nrlmsise00 neutral atmosphere model.

! the arguments are defined as follows:

!   iyr    the julian year
!   iday   the day of the year
!   hr     the local time in hours
!   glong  the geocentric longitude in degrees east
!   iyd    the year and day in the form yydd
!   secut  the universal time in seconds

    iyd    = 1000 * mod(iyr,100) + iday
    hrut   = hr - glong /15.

    do while ( hrut < 0.  )
        hrut = hrut + 24.
    enddo

    do while ( hrut >= 24. )
        hrut = hrut - 24.
    enddo

    secut  = hrut * 3600.

    return
    end subroutine msistim

! ***************************
! ***************************

!     get_gitm_data_i

! ***************************
! ***************************

    subroutine get_gitm_data_i(hrut,geolat,geolon,galt,un1,vn1 &
!    ,ttn1,tn21,to21,to11,tno1,tn4s1)
    ,ttn1,tn21,to21,to11,tno1,tn4s1,thyd1)

! IN: hrut,geolat,geolon,galt
! OUT: un1,vn1 are un(galt,geolat,geolon,hrut) and vn(...)
! a call to readunvn has previously loaded the /UNVN/

    use gitm_mod
    use grid_mod
    use variable_mod
!    use message_passing_mod
!    use parameter_mod

    character(len=80) :: str,chartime,file_name

! now to GITM lookup
! assume uniform angular grids
 
    dellat = lat(2)-lat(1)
    dellon = lon(2)-lon(1)
    fixlat = lat(nLats)
    fixlon = lon(1)

    slat = geolat
    slon = geolon
    salt = galt

    nalt = nAlts
    nlon = nLons
    nlat = nLats

! get indices in table for lat and long interpolation

      j = 1

      do while ( slat .ge. lat(j) .and. j .lt. nLats )
        ilat = j
        dellat = lat(j+1) - lat(j)
        j = j + 1
      enddo

      if ( slon < 0. ) slon = slon + 360.

      slon0 = slon - fixlon ! + 180.
      ilon = ifix(slon0/dellon) + 1


      if ( slon .le. lon(1) .and. slon .ge. 0 ) then
        slon = slon + 360.
        ilon = nLons
      endif
        
      fx = (slat - lat(ilat))/dellat
      fy = (slon - lon(ilon))/dellon
      xf = 1.-fx
      yf = 1.-fy

!        print *,ilon,slon,glont(ilon),fy,yf

!      if (fx .gt. 1 .or. fx .lt. 0) &
!          print *,ilat,slat,lat(ilat),fx,xf

    if ( fx .gt. 1 ) fx = 1.
    if ( fx .lt. 0 ) fx = 0.

    xf = 1.-fx
    yf = 1.-fy

      ilonp1 = ilon+1
      if( ilonp1 .eq. nlon+1) ilonp1 = 1
      ilatp1 = ilat+1

! now interpolate in altitude

    kalt = nalt
    do k=1,nalt
        if( altitudei(ilon,ilat,k) >= salt ) then
            kalt = k-1
            goto 996
        endif
    enddo
    996 continue

    kaltp1 = kalt+1

!!$    if( (ilon < 1) .OR. (ilat > nlon) ) &
!!$    write(6,33) taskid,ilon,ilat,kalt,slon,slat,salt
!!$    if( (ilat < 1) .OR. (ilat >= nlat) ) &
!!$    write(6,33) taskid,ilon,ilat,kalt,slon,slat,salt
!!$    if( (kalt < 1)  ) &
!!$    write(6,33) taskid,ilon,ilat,kalt,slon,slat,salt
!!$    33 format('UNVN1: ',4i4,1p3e10.2)

    fz=-1.

    fz=-1.

    if (kaltp1 <= nalt ) then
        fz = (salt - altitudei(ilon,ilat,kalt))/( &
        altitudei(ilon,ilat,kaltp1)-altitudei(ilon,ilat,kalt))

    ! low
        un1lo = uni(ilon,ilat,kalt)*xf*yf + &
        uni(ilon,ilatp1,kalt)*fx*yf + &
        uni(ilonp1,ilat,kalt)*xf*fy + &
        uni(ilonp1,ilatp1,kalt)*fx*fy
        vn1lo = vni(ilon,ilat,kalt)*xf*yf + &
        vni(ilon,ilatp1,kalt)*fx*yf + &
        vni(ilonp1,ilat,kalt)*xf*fy + &
        vni(ilonp1,ilatp1,kalt)*fx*fy
        tn1lo = temperaturei(ilon,ilat,kalt)*xf*yf + &
        temperaturei(ilon,ilatp1,kalt)*fx*yf + &
        temperaturei(ilonp1,ilat,kalt)*xf*fy + &
        temperaturei(ilonp1,ilatp1,kalt)*fx*fy
        tn21lo = n2i(ilon,ilat,kalt)*xf*yf + &
        n2i(ilon,ilatp1,kalt)*fx*yf + &
        n2i(ilonp1,ilat,kalt)*xf*fy + &
        n2i(ilonp1,ilatp1,kalt)*fx*fy
        to21lo = o2i(ilon,ilat,kalt)*xf*yf + &
        o2i(ilon,ilatp1,kalt)*fx*yf + &
        o2i(ilonp1,ilat,kalt)*xf*fy + &
        o2i(ilonp1,ilatp1,kalt)*fx*fy
        to11lo = oi(ilon,ilat,kalt)*xf*yf + &
        oi(ilon,ilatp1,kalt)*fx*yf + &
        oi(ilonp1,ilat,kalt)*xf*fy + &
        oi(ilonp1,ilatp1,kalt)*fx*fy
        tno1lo = noi(ilon,ilat,kalt)*xf*yf + &
        noi(ilon,ilatp1,kalt)*fx*yf + &
        noi(ilonp1,ilat,kalt)*xf*fy + &
        noi(ilonp1,ilatp1,kalt)*fx*fy
        tn4s1lo = n4si(ilon,ilat,kalt)*xf*yf + &
        n4si(ilon,ilatp1,kalt)*fx*yf + &
        n4si(ilonp1,ilat,kalt)*xf*fy + &
        n4si(ilonp1,ilatp1,kalt)*fx*fy
        thyd1lo = hi(ilon,ilat,kalt)*xf*yf + &
        hi(ilon,ilatp1,kalt)*fx*yf + &
        hi(ilonp1,ilat,kalt)*xf*fy + &
        hi(ilonp1,ilatp1,kalt)*fx*fy
    ! high
        un1hi = uni(ilon,ilat,kaltp1)*xf*yf + &
        uni(ilon,ilatp1,kaltp1)*fx*yf + &
        uni(ilonp1,ilat,kaltp1)*xf*fy + &
        uni(ilonp1,ilatp1,kaltp1)*fx*fy
        vn1hi = vni(ilon,ilat,kaltp1)*xf*yf + &
        vni(ilon,ilatp1,kaltp1)*fx*yf + &
        vni(ilonp1,ilat,kaltp1)*xf*fy + &
        vni(ilonp1,ilatp1,kaltp1)*fx*fy
        tn1hi = temperaturei(ilon,ilat,kaltp1)*xf*yf + &
        temperaturei(ilon,ilatp1,kaltp1)*fx*yf + &
        temperaturei(ilonp1,ilat,kaltp1)*xf*fy + &
        temperaturei(ilonp1,ilatp1,kaltp1)*fx*fy
        tn21hi = n2i(ilon,ilat,kaltp1)*xf*yf + &
        n2i(ilon,ilatp1,kaltp1)*fx*yf + &
        n2i(ilonp1,ilat,kaltp1)*xf*fy + &
        n2i(ilonp1,ilatp1,kaltp1)*fx*fy
        to21hi = o2i(ilon,ilat,kaltp1)*xf*yf + &
        o2i(ilon,ilatp1,kaltp1)*fx*yf + &
        o2i(ilonp1,ilat,kaltp1)*xf*fy + &
        o2i(ilonp1,ilatp1,kaltp1)*fx*fy
        to11hi = oi(ilon,ilat,kaltp1)*xf*yf + &
        oi(ilon,ilatp1,kaltp1)*fx*yf + &
        oi(ilonp1,ilat,kaltp1)*xf*fy + &
        oi(ilonp1,ilatp1,kaltp1)*fx*fy
        tno1hi = noi(ilon,ilat,kaltp1)*xf*yf + &
        noi(ilon,ilatp1,kaltp1)*fx*yf + &
        noi(ilonp1,ilat,kaltp1)*xf*fy + &
        noi(ilonp1,ilatp1,kaltp1)*fx*fy
        tn4s1hi = n4si(ilon,ilat,kaltp1)*xf*yf + &
        n4si(ilon,ilatp1,kaltp1)*fx*yf + &
        n4si(ilonp1,ilat,kaltp1)*xf*fy + &
        n4si(ilonp1,ilatp1,kaltp1)*fx*fy
        thyd1hi = hi(ilon,ilat,kaltp1)*xf*yf + &
        hi(ilon,ilatp1,kaltp1)*fx*yf + &
        hi(ilonp1,ilat,kaltp1)*xf*fy + &
        hi(ilonp1,ilatp1,kaltp1)*fx*fy

        un1 = un1lo*(1.-fz) + un1hi*fz
        vn1 = vn1lo*(1.-fz) + vn1hi*fz
        ttn1 = tn1lo*(1.-fz) + tn1hi*fz
        tn21 = tn21lo*(1.-fz) + tn21hi*fz
        to21 = to21lo*(1.-fz) + to21hi*fz
        to11 = to11lo*(1.-fz) + to11hi*fz
        tno1 = tno1lo*(1.-fz) + tno1hi*fz
        tn4s1 = tn4s1lo*(1.-fz) + tn4s1hi*fz
        thyd1 = thyd1lo*(1.-fz) + thyd1hi*fz


    else
    ! above the GITM grid:
    ! un,vn,tn are constant in altitude
        un1lo = uni(ilon,ilat,kalt)*xf*yf + &
        uni(ilon,ilatp1,kalt)*fx*yf + &
        uni(ilonp1,ilat,kalt)*xf*fy + &
        uni(ilonp1,ilatp1,kalt)*fx*fy
        vn1lo = vni(ilon,ilat,kalt)*xf*yf + &
        vni(ilon,ilatp1,kalt)*fx*yf + &
        vni(ilonp1,ilat,kalt)*xf*fy + &
        vni(ilonp1,ilatp1,kalt)*fx*fy
        tn1lo = temperaturei(ilon,ilat,kalt)*xf*yf + &
        temperaturei(ilon,ilatp1,kalt)*fx*yf + &
        temperaturei(ilonp1,ilat,kalt)*xf*fy + &
        temperaturei(ilonp1,ilatp1,kalt)*fx*fy
        un1 = un1lo
        vn1 = vn1lo
        ttn1 = tn1lo

    ! n2,o2,o,no,n4s,h are extrapolated exponentially

        tn21lo = n2i(ilon,ilat,kalt)*xf*yf + &
        n2i(ilon,ilatp1,kalt)*fx*yf + &
        n2i(ilonp1,ilat,kalt)*xf*fy + &
        n2i(ilonp1,ilatp1,kalt)*fx*fy
        to21lo = o2i(ilon,ilat,kalt)*xf*yf + &
        o2i(ilon,ilatp1,kalt)*fx*yf + &
        o2i(ilonp1,ilat,kalt)*xf*fy + &
        o2i(ilonp1,ilatp1,kalt)*fx*fy
        to11lo = oi(ilon,ilat,kalt)*xf*yf + &
        oi(ilon,ilatp1,kalt)*fx*yf + &
        oi(ilonp1,ilat,kalt)*xf*fy + &
        oi(ilonp1,ilatp1,kalt)*fx*fy
        tno1lo = noi(ilon,ilat,kalt)*xf*yf + &
        noi(ilon,ilatp1,kalt)*fx*yf + &
        noi(ilonp1,ilat,kalt)*xf*fy + &
        noi(ilonp1,ilatp1,kalt)*fx*fy
        tn4s1lo = n4si(ilon,ilat,kalt)*xf*yf + &
        n4si(ilon,ilatp1,kalt)*fx*yf + &
        n4si(ilonp1,ilat,kalt)*xf*fy + &
        n4si(ilonp1,ilatp1,kalt)*fx*fy
        thyd1lo = hi(ilon,ilat,kalt)*xf*yf + &
        hi(ilon,ilatp1,kalt)*fx*yf + &
        hi(ilonp1,ilat,kalt)*xf*fy + &
        hi(ilonp1,ilatp1,kalt)*fx*fy

    ! use highest grid point (kalt) and 5 from top
        kaltm = kalt-5
    ! note that though labelled hi, the following are at lower alt

        tn21hi = n2i(ilon,ilat,kaltm)*xf*yf + &
        n2i(ilon,ilatp1,kaltm)*fx*yf + &
        n2i(ilonp1,ilat,kaltm)*xf*fy + &
        n2i(ilonp1,ilatp1,kaltm)*fx*fy
        to21hi = o2i(ilon,ilat,kaltm)*xf*yf + &
        o2i(ilon,ilatp1,kaltm)*fx*yf + &
        o2i(ilonp1,ilat,kaltm)*xf*fy + &
        o2i(ilonp1,ilatp1,kaltm)*fx*fy
        to11hi = oi(ilon,ilat,kaltm)*xf*yf + &
        oi(ilon,ilatp1,kaltm)*fx*yf + &
        oi(ilonp1,ilat,kaltm)*xf*fy + &
        oi(ilonp1,ilatp1,kaltm)*fx*fy
        tno1hi = noi(ilon,ilat,kaltm)*xf*yf + &
        noi(ilon,ilatp1,kaltm)*fx*yf + &
        noi(ilonp1,ilat,kaltm)*xf*fy + &
        noi(ilonp1,ilatp1,kaltm)*fx*fy
        tn4s1hi = n4si(ilon,ilat,kaltm)*xf*yf + &
        n4si(ilon,ilatp1,kaltm)*fx*yf + &
        n4si(ilonp1,ilat,kaltm)*xf*fy + &
        n4si(ilonp1,ilatp1,kaltm)*fx*fy
        thyd1hi = hi(ilon,ilat,kaltm)*xf*yf + &
        hi(ilon,ilatp1,kaltm)*fx*yf + &
        hi(ilonp1,ilat,kaltm)*xf*fy + &
        hi(ilonp1,ilatp1,kaltm)*fx*fy

        hite2 = altitudei(ilon,ilat,kalt)
        hite1 = altitudei(ilon,ilat,kaltm)
               
        y2_den = alog( tn21lo )
        y1_den = alog( tn21hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        tn21 = exp(z_den)
               
        y2_den = alog( to21lo )
        y1_den = alog( to21hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        to21 = exp(z_den)

        y2_den = alog( to11lo )
        y1_den = alog( to11hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        to11 = exp(z_den)
               
        y2_den = alog( tno1lo )
        y1_den = alog( tno1hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        tno1 = exp(z_den)
               
        y2_den = alog( tn4s1lo )
        y1_den = alog( tn4s1hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        tn4s1 = exp(z_den)
               
        y2_den = alog( thyd1lo )
        y1_den = alog( thyd1hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        thyd1 = exp(z_den)

    endif

!!$    if (isnan(tn2int)) then
!!$    print *,'in gitm_get',&
!!$    hrut,geolat,geolon,galt,untint,vntint &
!!$    ,ttnint,tn2int,to2int,to1int,tnoint,tn4sint,theint
!!$    print *,'in gitm_get',&
!!$    hrut,tn21,to22
!!$    endif



! these are in m/sec

    end subroutine get_gitm_data_i


! ***************************
! ***************************

!     get_gitm_data_f

! ***************************
! ***************************

    subroutine get_gitm_data_f(hrut,geolat,geolon,galt,un1,vn1 &
!    ,ttn1,tn21,to21,to11,tno1,tn4s1)
    ,ttn1,tn21,to21,to11,tno1,tn4s1,thyd1)

! IN: hrut,geolat,geolon,galt
! OUT: un1,vn1 are un(galt,geolat,geolon,hrut) and vn(...)
! a call to readunvn has previously loaded the /UNVN/

    use gitm_mod
    use grid_mod
    use variable_mod
    use message_passing_mod
    use parameter_mod

    character(len=80) :: str,chartime,file_name

! now to GITM lookup
! assume uniform angular grids
 
    dellat = lat(2)-lat(1)
    dellon = lon(2)-lon(1)
    fixlat = lat(nLats)
    fixlon = lon(1)

    slat = geolat
    slon = geolon
    salt = galt

    nalt = nAlts
    nlon = nLons
    nlat = nLats

! get indices in table for lat and long interpolation

      j = 1

      do while ( slat .ge. lat(j) .and. j .lt. nLats )
        ilat = j
        dellat = lat(j+1) - lat(j)
        j = j + 1
      enddo

      if ( slon < 0. ) slon = slon + 360.

      slon0 = slon - fixlon ! + 180.
      ilon = ifix(slon0/dellon) + 1

      if ( slon .le. lon(1) .and. slon .ge. 0 ) then
        slon = slon + 360.
        ilon = nLons
      endif

!     if ( taskid == 1 ) print *,'i:',slon,slon0,dellon,ilon
!     stop
        
      fx = (slat - lat(ilat))/dellat
      fy = (slon - lon(ilon))/dellon
      xf = 1.-fx
      yf = 1.-fy

!        print *,ilon,slon,glont(ilon),fy,yf

!      if (fx .gt. 1 .or. fx .lt. 0) &
!          print *,ilat,slat,lat(ilat),fx,xf

    if ( fx .gt. 1 ) fx = 1.
    if ( fx .lt. 0 ) fx = 0.

    xf = 1.-fx
    yf = 1.-fy

      ilonp1 = ilon+1
      if( ilonp1 .eq. nlon+1) ilonp1 = 1
      ilatp1 = ilat+1

! now interpolate in altitude

    kalt = nalt
    do k=1,nalt
        if( altitudef(ilon,ilat,k) >= salt ) then
            kalt = k-1
            goto 996
        endif
    enddo
    996 continue

    kaltp1 = kalt+1

!!$    if( (ilon < 1) .OR. (ilat > nlon) ) &
!!$    write(6,33) taskid,ilon,ilat,kalt,slon,slat,salt
!!$    if( (ilat < 1) .OR. (ilat >= nlat) ) &
!!$    write(6,33) taskid,ilon,ilat,kalt,slon,slat,salt
!!$    if( (kalt < 1)  ) &
!!$    write(6,33) taskid,ilon,ilat,kalt,slon,slat,salt
!!$    33 format('UNVN1: ',4i4,1p3e10.2)

    fz=-1.

    if (kaltp1 <= nalt ) then
        fz = (salt - altitudef(ilon,ilat,kalt))/( &
        altitudef(ilon,ilat,kaltp1)-altitudef(ilon,ilat,kalt))

    ! low
        un1lo = unf(ilon,ilat,kalt)*xf*yf + &
        unf(ilon,ilatp1,kalt)*fx*yf + &
        unf(ilonp1,ilat,kalt)*xf*fy + &
        unf(ilonp1,ilatp1,kalt)*fx*fy
        vn1lo = vnf(ilon,ilat,kalt)*xf*yf + &
        vnf(ilon,ilatp1,kalt)*fx*yf + &
        vnf(ilonp1,ilat,kalt)*xf*fy + &
        vnf(ilonp1,ilatp1,kalt)*fx*fy
        tn1lo = temperaturef(ilon,ilat,kalt)*xf*yf + &
        temperaturef(ilon,ilatp1,kalt)*fx*yf + &
        temperaturef(ilonp1,ilat,kalt)*xf*fy + &
        temperaturef(ilonp1,ilatp1,kalt)*fx*fy
        tn21lo = n2f(ilon,ilat,kalt)*xf*yf + &
        n2f(ilon,ilatp1,kalt)*fx*yf + &
        n2f(ilonp1,ilat,kalt)*xf*fy + &
        n2f(ilonp1,ilatp1,kalt)*fx*fy
        to21lo = o2f(ilon,ilat,kalt)*xf*yf + &
        o2f(ilon,ilatp1,kalt)*fx*yf + &
        o2f(ilonp1,ilat,kalt)*xf*fy + &
        o2f(ilonp1,ilatp1,kalt)*fx*fy
        to11lo = of(ilon,ilat,kalt)*xf*yf + &
        of(ilon,ilatp1,kalt)*fx*yf + &
        of(ilonp1,ilat,kalt)*xf*fy + &
        of(ilonp1,ilatp1,kalt)*fx*fy
        tno1lo = nof(ilon,ilat,kalt)*xf*yf + &
        nof(ilon,ilatp1,kalt)*fx*yf + &
        nof(ilonp1,ilat,kalt)*xf*fy + &
        nof(ilonp1,ilatp1,kalt)*fx*fy
        tn4s1lo = n4sf(ilon,ilat,kalt)*xf*yf + &
        n4sf(ilon,ilatp1,kalt)*fx*yf + &
        n4sf(ilonp1,ilat,kalt)*xf*fy + &
        n4sf(ilonp1,ilatp1,kalt)*fx*fy
        thyd1lo = hf(ilon,ilat,kalt)*xf*yf + &
        hf(ilon,ilatp1,kalt)*fx*yf + &
        hf(ilonp1,ilat,kalt)*xf*fy + &
        hf(ilonp1,ilatp1,kalt)*fx*fy
    ! high
        un1hi = unf(ilon,ilat,kaltp1)*xf*yf + &
        unf(ilon,ilatp1,kaltp1)*fx*yf + &
        unf(ilonp1,ilat,kaltp1)*xf*fy + &
        unf(ilonp1,ilatp1,kaltp1)*fx*fy
        vn1hi = vnf(ilon,ilat,kaltp1)*xf*yf + &
        vnf(ilon,ilatp1,kaltp1)*fx*yf + &
        vnf(ilonp1,ilat,kaltp1)*xf*fy + &
        vnf(ilonp1,ilatp1,kaltp1)*fx*fy
        tn1hi = temperaturef(ilon,ilat,kaltp1)*xf*yf + &
        temperaturef(ilon,ilatp1,kaltp1)*fx*yf + &
        temperaturef(ilonp1,ilat,kaltp1)*xf*fy + &
        temperaturef(ilonp1,ilatp1,kaltp1)*fx*fy
        tn21hi = n2f(ilon,ilat,kaltp1)*xf*yf + &
        n2f(ilon,ilatp1,kaltp1)*fx*yf + &
        n2f(ilonp1,ilat,kaltp1)*xf*fy + &
        n2f(ilonp1,ilatp1,kaltp1)*fx*fy
        to21hi = o2f(ilon,ilat,kaltp1)*xf*yf + &
        o2f(ilon,ilatp1,kaltp1)*fx*yf + &
        o2f(ilonp1,ilat,kaltp1)*xf*fy + &
        o2f(ilonp1,ilatp1,kaltp1)*fx*fy
        to11hi = of(ilon,ilat,kaltp1)*xf*yf + &
        of(ilon,ilatp1,kaltp1)*fx*yf + &
        of(ilonp1,ilat,kaltp1)*xf*fy + &
        of(ilonp1,ilatp1,kaltp1)*fx*fy
        tno1hi = nof(ilon,ilat,kaltp1)*xf*yf + &
        nof(ilon,ilatp1,kaltp1)*fx*yf + &
        nof(ilonp1,ilat,kaltp1)*xf*fy + &
        nof(ilonp1,ilatp1,kaltp1)*fx*fy
        tn4s1hi = n4sf(ilon,ilat,kaltp1)*xf*yf + &
        n4sf(ilon,ilatp1,kaltp1)*fx*yf + &
        n4sf(ilonp1,ilat,kaltp1)*xf*fy + &
        n4sf(ilonp1,ilatp1,kaltp1)*fx*fy
        thyd1hi = hf(ilon,ilat,kaltp1)*xf*yf + &
        hf(ilon,ilatp1,kaltp1)*fx*yf + &
        hf(ilonp1,ilat,kaltp1)*xf*fy + &
        hf(ilonp1,ilatp1,kaltp1)*fx*fy

        un1 = un1lo*(1.-fz) + un1hi*fz
        vn1 = vn1lo*(1.-fz) + vn1hi*fz
        ttn1 = tn1lo*(1.-fz) + tn1hi*fz
        tn21 = tn21lo*(1.-fz) + tn21hi*fz
        to21 = to21lo*(1.-fz) + to21hi*fz
        to11 = to11lo*(1.-fz) + to11hi*fz
        tno1 = tno1lo*(1.-fz) + tno1hi*fz
        tn4s1 = tn4s1lo*(1.-fz) + tn4s1hi*fz
        thyd1 = thyd1lo*(1.-fz) + thyd1hi*fz


    else
    ! above the GITM grid:
    ! un,vn,tn are constant in altitude
        un1lo = unf(ilon,ilat,kalt)*xf*yf + &
        unf(ilon,ilatp1,kalt)*fx*yf + &
        unf(ilonp1,ilat,kalt)*xf*fy + &
        unf(ilonp1,ilatp1,kalt)*fx*fy
        vn1lo = vnf(ilon,ilat,kalt)*xf*yf + &
        vnf(ilon,ilatp1,kalt)*fx*yf + &
        vnf(ilonp1,ilat,kalt)*xf*fy + &
        vnf(ilonp1,ilatp1,kalt)*fx*fy
        tn1lo = temperaturef(ilon,ilat,kalt)*xf*yf + &
        temperaturef(ilon,ilatp1,kalt)*fx*yf + &
        temperaturef(ilonp1,ilat,kalt)*xf*fy + &
        temperaturef(ilonp1,ilatp1,kalt)*fx*fy
        un1 = un1lo
        vn1 = vn1lo
        ttn1 = tn1lo

    ! n2,o2,o,no,n4s,h are extrapolated exponentially

        tn21lo = n2f(ilon,ilat,kalt)*xf*yf + &
        n2f(ilon,ilatp1,kalt)*fx*yf + &
        n2f(ilonp1,ilat,kalt)*xf*fy + &
        n2f(ilonp1,ilatp1,kalt)*fx*fy
        to21lo = o2f(ilon,ilat,kalt)*xf*yf + &
        o2f(ilon,ilatp1,kalt)*fx*yf + &
        o2f(ilonp1,ilat,kalt)*xf*fy + &
        o2f(ilonp1,ilatp1,kalt)*fx*fy
        to11lo = of(ilon,ilat,kalt)*xf*yf + &
        of(ilon,ilatp1,kalt)*fx*yf + &
        of(ilonp1,ilat,kalt)*xf*fy + &
        of(ilonp1,ilatp1,kalt)*fx*fy
        tno1lo = nof(ilon,ilat,kalt)*xf*yf + &
        nof(ilon,ilatp1,kalt)*fx*yf + &
        nof(ilonp1,ilat,kalt)*xf*fy + &
        nof(ilonp1,ilatp1,kalt)*fx*fy
        tn4s1lo = n4sf(ilon,ilat,kalt)*xf*yf + &
        n4sf(ilon,ilatp1,kalt)*fx*yf + &
        n4sf(ilonp1,ilat,kalt)*xf*fy + &
        n4sf(ilonp1,ilatp1,kalt)*fx*fy
        thyd1lo = hf(ilon,ilat,kalt)*xf*yf + &
        hf(ilon,ilatp1,kalt)*fx*yf + &
        hf(ilonp1,ilat,kalt)*xf*fy + &
        hf(ilonp1,ilatp1,kalt)*fx*fy

    ! use highest grid point (kalt) and 5 from top
        kaltm = kalt-5
    ! note that though labelled hi, the following are at lower alt

        tn21hi = n2f(ilon,ilat,kaltm)*xf*yf + &
        n2f(ilon,ilatp1,kaltm)*fx*yf + &
        n2f(ilonp1,ilat,kaltm)*xf*fy + &
        n2f(ilonp1,ilatp1,kaltm)*fx*fy
        to21hi = o2f(ilon,ilat,kaltm)*xf*yf + &
        o2f(ilon,ilatp1,kaltm)*fx*yf + &
        o2f(ilonp1,ilat,kaltm)*xf*fy + &
        o2f(ilonp1,ilatp1,kaltm)*fx*fy
        to11hi = of(ilon,ilat,kaltm)*xf*yf + &
        of(ilon,ilatp1,kaltm)*fx*yf + &
        of(ilonp1,ilat,kaltm)*xf*fy + &
        of(ilonp1,ilatp1,kaltm)*fx*fy
        tno1hi = nof(ilon,ilat,kaltm)*xf*yf + &
        nof(ilon,ilatp1,kaltm)*fx*yf + &
        nof(ilonp1,ilat,kaltm)*xf*fy + &
        nof(ilonp1,ilatp1,kaltm)*fx*fy
        tn4s1hi = n4sf(ilon,ilat,kaltm)*xf*yf + &
        n4sf(ilon,ilatp1,kaltm)*fx*yf + &
        n4sf(ilonp1,ilat,kaltm)*xf*fy + &
        n4sf(ilonp1,ilatp1,kaltm)*fx*fy
        thyd1hi = hf(ilon,ilat,kaltm)*xf*yf + &
        hf(ilon,ilatp1,kaltm)*fx*yf + &
        hf(ilonp1,ilat,kaltm)*xf*fy + &
        hf(ilonp1,ilatp1,kaltm)*fx*fy

        hite2 = altitudef(ilon,ilat,kalt)
        hite1 = altitudef(ilon,ilat,kaltm)
               
        y2_den = alog( tn21lo )
        y1_den = alog( tn21hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        tn21 = exp(z_den)
               
        y2_den = alog( to21lo )
        y1_den = alog( to21hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        to21 = exp(z_den)

        y2_den = alog( to11lo )
        y1_den = alog( to11hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        to11 = exp(z_den)
               
        y2_den = alog( tno1lo )
        y1_den = alog( tno1hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        tno1 = exp(z_den)
               
        y2_den = alog( tn4s1lo )
        y1_den = alog( tn4s1hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        tn4s1 = exp(z_den)
               
        y2_den = alog( thyd1lo )
        y1_den = alog( thyd1hi )
        slope = (y2_den-y1_den)/(hite2-hite1)
        z_den = y2_den + slope*(salt - hite2)
        thyd1 = exp(z_den)

    endif

!!$    if (isnan(tn2int)) then
!!$    print *,'in gitm_get',&
!!$    hrut,geolat,geolon,galt,untint,vntint &
!!$    ,ttnint,tn2int,to2int,to1int,tnoint,tn4sint,theint
!!$    print *,'in gitm_get',&
!!$    hrut,tn21,to22
!!$    endif

! these are in m/sec

    end subroutine get_gitm_data_f


  subroutine read_gitm_for_sami3(nfile)

    use gitm_mod
    use message_passing_mod

    character (len=80) :: FileNameIn,str,chartime,dir

    write(chartime,'(i5)') 10000+nfile
    str=chartime(2:)
    dir='./gitm_data/'
    FileNameIn  = trim(dir)//'gitm_'//trim(str)//'.inp'
    print *,'file_name =',FileNameIn 
    close(101)
    open(101,file=FileNameIn,form="unformatted")

  read(101) inLons, inLats, inAlts, inFiles
  read(101) lat
  read(101) lon
  read(101) altitudei
  read(101) uni
  read(101) vni
  read(101) temperaturei
  read(101) n2i
  read(101) o2i
  read(101) oi
  read(101) n4si
  read(101) noi
  read(101) hi

! convert to km

  altitudei = altitudei / 1.e3

! convert to cm^-3

  n2i       = n2i / 1.e6
  o2i       = o2i / 1.e6
  oi        = oi  / 1.e6
  n4si      = n4si/ 1.e6
  noi       = noi / 1.e6
  hi        = hi  / 1.e6

!!$  print *,'uni,vni ',taskid,uni(30,50,20),vni(30,50,20)
!!$  print *,'o,o2 ',taskid,oi(30,50,20),o2i(30,50,20)

    write(chartime,'(i5)') 10000+(nfile+1)
    str=chartime(2:)
    dir='./gitm_data/'
    FileNameIn  = trim(dir)//'gitm_'//trim(str)//'.inp'
    print *,'file_name =',FileNameIn 
    close(101)
    open(101,file=FileNameIn,form="unformatted")

  read(101) inLons, inLats, inAlts, inFiles
  read(101) lat
  read(101) lon
  read(101) altitudef
  read(101) unf
  read(101) vnf
  read(101) temperaturef
  read(101) n2f
  read(101) o2f
  read(101) of
  read(101) n4sf
  read(101) nof
  read(101) hf

! convert to km

  altitudef = altitudef / 1.e3

! convert to cm^-3

  n2f       = n2f / 1.e6
  o2f       = o2f / 1.e6
  of        = of  / 1.e6
  n4sf      = n4sf/ 1.e6
  nof       = nof / 1.e6
  hf        = hf  / 1.e6

  

end subroutine read_gitm_for_sami3


