
!******************************************
!******************************************

!            chemrate

!******************************************
!******************************************

!     chemical producation and loss rates
!     bb: bailley and balan (red book, 1996)

    subroutine chemrate ( chrate,nfl,nll )

    use parameter_mod
    use variable_mod

    real :: chrate(nz,nchem)

    do iz = 1,nz

        ti300o = ti(iz,nfl,nll,ptop) / 300.

    ! h+ + o --> o+ + h (bb)

        chrate (iz,1) = 2.2e-11 &
        * sqrt( ti(iz,nfl,nll,pthp) )

    ! he+ + n2 --> n2+ + he (bb)

        chrate (iz,2) = 3.5e-10

    ! he+ + n2 --> n+ + n + he (schunk)

        chrate (iz,3) = 8.5e-10

    ! he+ + o2 --> o+ + o + he (bb)

        chrate (iz,4) = 8.0e-10

    ! he+ + o2 --> o2+ + he

        chrate (iz,5) = 2.0e-10

    ! n+ + o2 --> no+ + o  (schunk)

        chrate (iz,6) = 2.0e-10

    ! n+ + o2 --> o2+ + n(2d) (schunk)

        chrate (iz,7) = 4.0e-10

    ! n+ + 0 --> o+ + n

        chrate (iz,8) = 1.0e-12

    ! n+ + no --> no+ + o (schunk)

        chrate (iz,9) = 2.0e-11

    ! o+ + h --> h+ + o   (bb)

        chrate(iz,10) = 2.5e-11 &
        * sqrt( tn(iz,nfl,nll) )
         
    ! o+ + n2 --> no+ + n (bb)

        chrate(iz,11) = 1.533e-12 - &
        &              5.920e-13 * ti300o + &
        &              8.600e-14 * ti300o ** 2

        if ( ti(iz,nfl,nll,ptop) > 1700 ) &
        chrate(iz,11) = 2.730e-12 - &
        &                     1.155e-12 * ti300o + &
        &                     1.483e-13 * ti300o ** 2

    ! o+ + o2 --> o2+ + o

        chrate(iz,12) = 2.820e-11 - &
        &              7.740e-12 * ti300o + &
        &              1.073e-12 * ti300o ** 2 - &
        &              5.170e-14 * ti300o ** 3 + &
        &              9.650e-16 * ti300o ** 4

    ! o+ + no --> no+ + o

        chrate(iz,13) = 1.0e-12

    ! n2+ + o --> no+ + n(2d) (bb)

        chrate(iz,14) = 1.4e-10 / ti300o ** .44

    ! n2+ + o2 --> o2+ + n2 (schunk)

        chrate(iz,15) = 5.0e-11 / sqrt( ti300o )

    ! n2+ + o2 --> no+ + no

        chrate(iz,16) = 1.0e-14

    ! n2+ + no --> no+ + n2 (schunk)

        chrate(iz,17) = 3.3e-10

    ! o2+ + n --> no+ + o (schunk)

        chrate(iz,18) = 1.2e-10

    ! o2+ + n(2d) --> n+ + o2

        chrate(iz,19) = 2.5e-10

    ! o2+ + no --> no+ + o2 (bb)

        chrate(iz,20) = 4.4e-10

    ! o2+ + n2 --> no+ + no (schunk)

        chrate(iz,21) = 5.0e-16

    enddo

    return
    end subroutine chemrate

!******************************************
!******************************************

!            recorate

!******************************************
!******************************************

! recombination rates
! bb: bailley and balan (red book, 1996)

    subroutine recorate ( relossr,nfl,nll )

    use parameter_mod
    use variable_mod

    real :: relossr(nz,nion)

    do iz = 1,nz

        te300 = te(iz,nfl,nll) / 300.

        relossr(iz,pthp)  = 4.43e-12 / te300 ** .7
        relossr(iz,pthep) = relossr(iz,pthp)
        relossr(iz,ptnp)  = relossr(iz,pthp)
        relossr(iz,ptop)  = relossr(iz,pthp)
        relossr(iz,ptn2p) = 1.8e-7 / te300 ** 0.39     !   (schunk)
        relossr(iz,ptnop) = 4.2e-7 / te300 ** 0.85     !   (bb)
        relossr(iz,pto2p) = 1.6e-7 / te300 ** 0.55     !   (schunk)

    enddo

    return
    end subroutine recorate


!******************************************
!******************************************

!          chempl

!******************************************
!******************************************
            
! chemical loss (chloss) and production (chprod)

! chrate: chemical reaction rates calculated in chemrate
! ichem: input data file showing loss, neutral, production
!        species for each reaction

    subroutine chempl ( chrate,chloss,chprod,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use chemistry_mod

    real :: chrate(nz,nchem),chloss(nz,nion),chprod(nz,nion)

    do i = nion1,nion2
        do iz = 1,nz
            chloss(iz,i)   = 0.
            chprod(iz,i)   = 0.
        enddo
    enddo

    do k = 1,nchem
        il = ichem(k,1) ! ion species (reacting) loss
        in = ichem(k,2) ! neutral species reacting
        ip = ichem(k,3) ! ion species produced
        do iz = 1,nz
            chem  = denn(iz,nfl,nll,in) * chrate(iz,k)
            tdeni = deni(iz,nfl,nll,il) * chem
            chloss(iz,il) = tdeni + chloss(iz,il)
            chprod(iz,ip) = tdeni + chprod(iz,ip)
        enddo
    enddo

    return
    end subroutine chempl

