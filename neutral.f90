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

    do j = 1,nf
        do i = 1,nz
            glonsij = glons(i,j,nll)
            if ( lcr ) glonsij = glons0(i,j,nll) - hruti * 15.
!            if ( lcr ) glonsij = glons0(i,j,nll) - (hruti - hrinit) * 15.
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
       if (taskid.eq.1 .and. i.eq.nz/2 .and. j.eq.20 &
                      .and. nll.eq.nl/2 ) &
          print *,'d = ',d
        enddo
    enddo
    do j = 1,nf
        do i = 1,nz
            app(1)   = ap
            app(2)   = ap
            glonsij = glons(i,j,nll)
            if ( lcr ) glonsij = glons0(i,j,nll) - hruti * 15.
!            if ( lcr ) glonsij = glons0(i,j,nll) - (hruti - hrinit) * 15.
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
        enddo
    enddo

    hrutf  = hr + .25

    do j = 1,nf
        do i = 1,nz
            glonsij = glons(i,j,nll)
            if ( lcr ) glonsij = glons0(i,j,nll) - hrutf * 15.
!            if ( lcr ) glonsij = glons0(i,j,nll) - (hrutf - hrinit) * 15.
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
            if ( lcr ) glonsij = glons0(i,j,nll) - hrutf * 15.
!            if ( lcr ) glonsij = glons0(i,j,nll) - (hrutf - hrinit) * 15.
            glonsij = mod(glonsij,360.)
            hrl   = mod(hrutf + glonsij / 15.,24.)
            call msistim ( int(year),int(day),hrl,glonsij,iyd,sec )
            if(lhwm93)  call gws5 ( iyd,sec,alts(i,j,nll), &
                                    glats(i,j,nll),glonsij, &
                                    hrl,fbar,f10p7,app,whm93        )
            if(lhwm14)  call hwm14 ( iyd,sec,alts(i,j,nll), &
                                     glats(i,j,nll),glonsij, &
                                     hrl,fbar,f10p7,app,whm93        )

            vf(i,j,nll)   = 100. * whm93(1) * tvn0
            uf(i,j,nll)   = 100. * whm93(2) * tvn0
            wf(i,j,nll)   = vw   * tvn0

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
    use misc_mod

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
                u1(i,j,n) = u(i,j,n)
                u2(i,j,n) = v(i,j,n)
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
