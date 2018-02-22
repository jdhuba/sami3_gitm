!******************************************
!******************************************

!             EXB

!******************************************
!******************************************

    subroutine exb(hrut,phi)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use exb_mod
    use misc_mod
    use grid_mod

    real :: denic(nz,nf,nl,nion)
    real :: tic(nz,nf,nl,nion)
    real :: tec(nz,nf,nl)
    real :: fluxnp(nz,nfp1,nl,nion),fluxtp(nz,nfp1,nl,nion)
    real :: fluxtep(nz,nfp1,nl)
    real :: fluxns(nz,nf,nl,nion),fluxts(nz,nf,nl,nion)
    real :: fluxtes(nz,nf,nl)
    real :: fluxnh(nz,nf,nl,nion),fluxth(nz,nf,nl,nion)
    real :: fluxteh(nz,nf,nl)
    real :: vexb_p(nzp1,nfp1,nlp1),vexb_h(nzp1,nfp1,nlp1)
    real :: param(2)
    real :: phi(nnx,nny)

! define the e x b drift

    param(1) = day
    param(2) = f10p7
    nzh      = nz / 2

! note: modification of vexb because of field line variation
!       uses cos^3/sqrt(1.+3sin^2) instead of
!       uses sin^3/sqrt(1.+3cos^2) because
!       blats = 0 at the magnetic equator
!       (as opposed to pi/2 in spherical coordinates)

    call vexb_phi(phi)
    call current_jp_jphi

! here we add in exb drift from the potential

! reduce e x b velocities below alt_crit
! with 20 km decay (dela)

! and kill above some high altitude

    alt_crit_high   =  pcrit * re
    dela_high     = 2. * re
    dela = 20.

!     vexbp

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                vexbp(i,j,k) = vexbp_phi(i,j,k)
                if ( baltp(i,j,k) < alt_crit ) then
                    arg0 = ( alt_crit - baltp(i,j,k) ) / dela
                    fac = exp(-arg0*arg0)
                    vexbp(i,j,k) = vexbp(i,j,k) * fac
                endif
                if ( baltp(i,j,k) > alt_crit_high ) then
                    arg0 = ( abs(alt_crit_high - baltp(i,j,k)) ) / dela_high
                    fac = exp(-arg0*arg0)
                    vexbp(i,j,k) = vexbp(i,j,k) * fac
                endif
            enddo
        enddo
    enddo

!     vexbs

    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1
                vexbs(i,j,k) = vexbs_phi(i,j,k)
                if ( baltp(i,j,k) < alt_crit ) then
                    arg0 = ( alt_crit - baltp(i,j,k) ) / dela
                    fac = exp(-arg0*arg0)
                    vexbs(i,j,k) = vexbs(i,j,k) * fac
                endif
                if ( baltp(i,j,k) > alt_crit_high ) then
                    arg0 = ( abs(alt_crit_high - baltp(i,j,k)) ) / dela_high
                    fac = exp(-arg0*arg0)
                    vexbs(i,j,k) = vexbs(i,j,k) * fac
                endif
            enddo
        enddo
    enddo

!     vexbh

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz
                vexbh(i,j,k) = vexbh_phi(i,j,k)
                if ( baltp(i,j,k) < alt_crit ) then
                    arg0 = ( alt_crit - baltp(i,j,k) ) / dela
                    fac = exp(-arg0*arg0)
                    vexbh(i,j,k) = vexbh(i,j,k) * fac
                endif
                if ( baltp(i,j,k) > alt_crit_high ) then
                    arg0 = ( abs(alt_crit_high - baltp(i,j,k)) ) / dela_high
                    fac = exp(-arg0*arg0)
                    vexbh(i,j,k) = vexbh(i,j,k) * fac
                endif
            enddo
        enddo
    enddo

! limit e x b velocities

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                if (vexbp(i,j,k) > 0.) &
                vexbp(i,j,k) = amin1(vexbp(i,j,k),vexb_max)
                if (vexbp(i,j,k) < 0.) &
                vexbp(i,j,k) = amax1(vexbp(i,j,k),-vexb_max)
            enddo
        enddo
    enddo


    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1
                if (vexbs(i,j,k) > 0.) &
                vexbs(i,j,k) = amin1(vexbs(i,j,k),vexb_max)
                if (vexbs(i,j,k) < 0.) &
                vexbs(i,j,k) = amax1(vexbs(i,j,k),-vexb_max)
            enddo
        enddo
    enddo

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz
                if (vexbh(i,j,k) > 0.) &
                vexbh(i,j,k) = amin1(vexbh(i,j,k),vexb_max)
                if (vexbh(i,j,k) < 0.) &
                vexbh(i,j,k) = amax1(vexbh(i,j,k),-vexb_max)
            enddo
        enddo
    enddo


! output e x b velocities

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                u1p(i,j,k) = vexbp(i,j,k)
                u2s(i,j,k) = vexbs(i,j,k)
                u3h(i,j,k) = vexbh(i,j,k)
            enddo
        enddo
    enddo

! calculate conserved particle number: denic
! and 'conserved' temperature: tic,tec

    do ni = nion1,nion2
        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    denic(i,j,k,ni) = deni(i,j,k,ni) * vol(i,j,k)
                    tic(i,j,k,ni)   = ti(i,j,k,ni) * vol(i,j,k)
                enddo
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                tec(i,j,k)   = te(i,j,k) * vol(i,j,k)
            enddo
        enddo
    enddo

! calculate flux in p-direction at interface
! NOTE: neutral flux condition at outer boundary (JH 11/29/07)
! altered to consider NS pole densities

    do ni = nion1,nion2
        do k = 1,nl
            do j = 2,nf
                do i = 1,nz
                    if ( vexbp(i,j,k) >= 0 ) then
                        fluxnp(i,j,k,ni) = deni(i,j-1,k,ni) * vexbp(i,j,k)
                        fluxtp(i,j,k,ni) = ti(i,j-1,k,ni)   * vexbp(i,j,k)
                    else
                        fluxnp(i,j,k,ni) = deni(i,j,k,ni) * vexbp(i,j,k)
                        fluxtp(i,j,k,ni) = ti(i,j,k,ni)   * vexbp(i,j,k)
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 2,nf
            do i = 1,nz
                if ( vexbp(i,j,k) >= 0 ) then
                    fluxtep(i,j,k) = te(i,j-1,k) * vexbp(i,j,k)
                else
                    fluxtep(i,j,k) = te(i,j,k)   * vexbp(i,j,k)
                endif
            enddo
        enddo
    enddo

!      flux at nfp1 (near magnetic north/south poles)

    do ni = nion1,nion2
        do k = 1,nl
            j = nfp1
            do i = 1,nz
                if ( vexbp(i,j,k) >= 0 ) then
                    fluxnp(i,j,k,ni) = deni(i,j-1,k,ni) * vexbp(i,j,k)
                    fluxtp(i,j,k,ni) = ti(i,j-1,k,ni)   * vexbp(i,j,k)
                else
                    fluxnp(i,j,k,ni) = deni_mnp(i,ni) * vexbp(i,j,k)
                    fluxtp(i,j,k,ni) = ti_mnp(i,ni)   * vexbp(i,j,k)
                endif
            enddo
        enddo
    enddo

    do k = 1,nl
        j = nfp1
        do i = 1,nz
            if ( vexbp(i,j,k) >= 0 ) then
                fluxtep(i,j,k) = te(i,j-1,k) * vexbp(i,j,k)
            else
                fluxtep(i,j,k) = te_mnp(i)   * vexbp(i,j,k)
            endif
        enddo
    enddo

! calculate flux in s-direction at interface

    do ni = nion1,nion2
        do k = 1,nl
            do j = 1,nf
                do i = 2,nz
                    if ( vexbs(i,j,k) >= 0 ) then
                        fluxns(i,j,k,ni) = deni(i-1,j,k,ni) * vexbs(i,j,k)
                        fluxts(i,j,k,ni) = ti(i-1,j,k,ni)   * vexbs(i,j,k)
                    else
                        fluxns(i,j,k,ni) = deni(i,j,k,ni) * vexbs(i,j,k)
                        fluxts(i,j,k,ni) = ti(i,j,k,ni)   * vexbs(i,j,k)
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nf
            do i = 2,nz
                if ( vexbs(i,j,k) >= 0 ) then
                    fluxtes(i,j,k) = te(i-1,j,k) * vexbs(i,j,k)
                else
                    fluxtes(i,j,k) = te(i,j,k)   * vexbs(i,j,k)
                endif
            enddo
        enddo
    enddo

! calculate flux in h-direction at interface (k > 1)

    do ni = nion1,nion2
        do k = 2,nl
            do j = 1,nf
                do i = 1,nz
                    if ( vexbh(i,j,k) >= 0 ) then
                        fluxnh(i,j,k,ni) = deni(i,j,k-1,ni) * vexbh(i,j,k)
                        fluxth(i,j,k,ni) = ti(i,j,k-1,ni)   * vexbh(i,j,k)
                    else
                        fluxnh(i,j,k,ni) = deni(i,j,k,ni) * vexbh(i,j,k)
                        fluxth(i,j,k,ni) = ti(i,j,k,ni)   * vexbh(i,j,k)
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 2,nl
        do j = 1,nf
            do i = 1,nz
                if ( vexbh(i,j,k) >= 0 ) then
                    fluxteh(i,j,k) = te(i,j,k-1) * vexbh(i,j,k)
                else
                    fluxteh(i,j,k) = te(i,j,k)   * vexbh(i,j,k)
                endif
            enddo
        enddo
    enddo

!      calculate flux in h-direction at interface (k = 1)
!      (invoke periodic boundary condition)

!  modify: extrapolate

 
    do ni = nion1,nion2
        do j = 1,nf
            do i = 1,nz
               if ( vexbh(i,j,1) >= 0 ) then
                 del_blon = (blonp(i,j,3)-blonp(i,j,1)) / &
                           (blonp(i,j,2)-blonp(i,j,1))
                 xdeni    = deni(i,j,3,ni) - &
                           del_blon * ( deni(i,j,3,ni) - deni(i,j,2,ni) )
                 xti      = ti(i,j,3,ni) -   &
                           del_blon * ( ti(i,j,3,ni) - ti(i,j,2,ni) )
                 fluxnh(i,j,1,ni) = xdeni * vexbh(i,j,1)
                 fluxth(i,j,1,ni) = xti   * vexbh(i,j,1)

!                    fluxnh(i,j,1,ni) = deni(i,j,2,ni) * vexbh(i,j,1) ! zero gradient
!                    fluxth(i,j,1,ni) = ti(i,j,2,ni)   * vexbh(i,j,1) ! zero gradient
!                    fluxnh(i,j,1,ni) = deni(i,j,nl,ni) * vexbh(i,j,1) ! periodic
!                    fluxth(i,j,1,ni) = ti(i,j,nl,ni)   * vexbh(i,j,1) ! periodic
                else
                    fluxnh(i,j,1,ni) = deni(i,j,1,ni) * vexbh(i,j,1)
                    fluxth(i,j,1,ni) = ti(i,j,1,ni)   * vexbh(i,j,1)
                endif
            enddo
        enddo
    enddo

    do j = 1,nf
        do i = 1,nz
            if ( vexbh(i,j,1) >= 0 ) then
              del_blon = (blonp(i,j,3)-blonp(i,j,1)) / &
                         (blonp(i,j,2)-blonp(i,j,1))
              xte      = te(i,j,3) - &
                         del_blon * ( te(i,j,3) - te(i,j,2) )
              fluxteh(i,j,1) = xte * vexbh(i,j,1)
!                fluxteh(i,j,1) = te(i,j,2) * vexbh(i,j,1)  ! zero gradient
!                fluxteh(i,j,1) = te(i,j,nl) * vexbh(i,j,1) ! periodic 
            else
                fluxteh(i,j,1) = te(i,j,1)  * vexbh(i,j,1)
            endif
        enddo
    enddo

! update total particle number and density
! and temperatures
! NOTE: the temperature update is an approximation
!       (probably better than no update but, strictly
!       speaking, not exactly correct)

    do ni = nion1,nion2
        do k = 2,nlm1
            do j = 2,nf
                do i = 2,nzm1
                    denic(i,j,k,ni) = denic(i,j,k,ni) &
                    + dt * ( areap(i,j,k)   * fluxnp(i,j,k,ni) - &
                    areap(i,j+1,k) * fluxnp(i,j+1,k,ni) ) &
                    + dt * ( areas(i,j,k)   * fluxns(i,j,k,ni) - &
                    areas(i+1,j,k) * fluxns(i+1,j,k,ni) ) &
                    + dt * ( areah(i,j,k)   * fluxnh(i,j,k,ni) - &
                    areah(i,j,k+1) * fluxnh(i,j,k+1,ni) )
                    deni(i,j,k,ni)  = denic(i,j,k,ni) / vol(i,j,k)

                ! brazen fix
                    deni(i,j,k,ni)  = amax1(deni(i,j,k,ni),denmin)

                    tic(i,j,k,ni) = tic(i,j,k,ni) &
                    + dt * ( areap(i,j,k)   * fluxtp(i,j,k,ni) - &
                    areap(i,j+1,k) * fluxtp(i,j+1,k,ni) ) &
                    + dt * ( areas(i,j,k)   * fluxts(i,j,k,ni) - &
                    areas(i+1,j,k) * fluxts(i+1,j,k,ni) ) &
                    + dt * ( areah(i,j,k)   * fluxth(i,j,k,ni) - &
                    areah(i,j,k+1) * fluxth(i,j,k+1,ni) )
                    ti(i,j,k,ni)  = tic(i,j,k,ni) / vol(i,j,k)
                ! brazen fix
                    ti(i,j,k,ni)  = amax1(ti(i,j,k,ni),200.)
                    if (isnan(ti(i,j,k,ni))) then
                        print *,'Ti fixed',i,j,k,ni
                        ti(i,j,k,ni) = 200.
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 2,nlm1
        do j = 2,nf
            do i = 2,nzm1
                tec(i,j,k) = tec(i,j,k) &
                + dt * ( areap(i,j,k)   * fluxtep(i,j,k) - &
                areap(i,j+1,k) * fluxtep(i,j+1,k) ) &
                + dt * ( areas(i,j,k)   * fluxtes(i,j,k) - &
                areas(i+1,j,k) * fluxtes(i+1,j,k) ) &
                + dt * ( areah(i,j,k)   * fluxteh(i,j,k) - &
                areah(i,j,k+1) * fluxteh(i,j,k+1) )
                te(i,j,k)  = tec(i,j,k) / vol(i,j,k)
            ! brazen fix
                te(i,j,k)  = amax1(te(i,j,k),200.)
                if (te(i,j,k) < 0.) print *,'i,j,k', &
                i,j,k,te(i,j,k)
            enddo
        enddo
    enddo

!      for k = nl

    do ni = nion1,nion2
        do j = 2,nf
            do i = 2,nzm1
                k                = nl
                deni0            = deni(i,j,k,ni)
                ti0              = ti(i,j,k,ni)

                denic(i,j,nl,ni) = denic(i,j,nl,ni) &
                + dt * ( areap(i,j,nl)   * fluxnp(i,j,nl,ni) - &
                areap(i,j+1,nl) * fluxnp(i,j+1,nl,ni) ) &
                + dt * ( areas(i,j,nl)   * fluxns(i,j,nl,ni) - &
                areas(i+1,j,nl) * fluxns(i+1,j,nl,ni) ) &
                + dt * ( areah(i,j,nl)   * fluxnh(i,j,nl,ni) - &
                areah(i,j,1)    * fluxnh(i,j,1,ni) )

                deni(i,j,nl,ni)  = denic(i,j,nl,ni) / vol(i,j,nl)

                if ( deni(i,j,nl,ni) <= 0. ) &
                deni(i,j,nl,ni) = deni0
                tic(i,j,nl,ni) = tic(i,j,nl,ni) &
                + dt * ( areap(i,j,nl)   * fluxtp(i,j,nl,ni) - &
                areap(i,j+1,nl) * fluxtp(i,j+1,nl,ni) ) &
                + dt * ( areas(i,j,nl)   * fluxts(i,j,nl,ni) - &
                areas(i+1,j,nl) * fluxts(i+1,j,nl,ni) ) &
                + dt * ( areah(i,j,nl)   * fluxth(i,j,nl,ni) - &
                areah(i,j,1)    * fluxth(i,j,1,ni) )
                ti(i,j,nl,ni)  = tic(i,j,nl,ni) / vol(i,j,nl)
                if ( ti(i,j,nl,ni) <= 0. ) &
                ti(i,j,nl,ni) = ti0
            enddo
        enddo
    enddo

    do j = 2,nf
        do i = 2,nzm1
            te0 = te(i,j,nl)
            tec(i,j,nl) = tec(i,j,nl) &
            + dt * ( areap(i,j,nl)   * fluxtep(i,j,nl) - &
            areap(i,j+1,nl) * fluxtep(i,j+1,nl) ) &
            + dt * ( areas(i,j,nl)   * fluxtes(i,j,nl) - &
            areas(i+1,j,nl) * fluxtes(i+1,j,nl) ) &
            + dt * ( areah(i,j,nl)   * fluxteh(i,j,nl) - &
            areah(i,j,1)    * fluxteh(i,j,1) )
            te(i,j,nl)  = tec(i,j,nl) / vol(i,j,nl)
            if ( te(i,j,nl) <= 0. ) &
            te(i,j,nl) = te0
        enddo
    enddo


! fill cells at j = 1  with j = 2 

    do ni = nion1,nion2
        do k = 2,nlm1
            do i = 2,nzm1
                deni(i,1,k,ni)  = deni(i,2,k,ni)
            !             deni(i,nf,k,ni) = deni(i,nfm1,k,ni)
                ti(i,1,k,ni)    = ti(i,2,k,ni)
            !             ti(i,nf,k,ni)   = ti(i,nfm1,k,ni)
            enddo
        enddo
    enddo

    do k = 2,nlm1
        do i = 2,nzm1
            te(i,1,k)    = te(i,2,k)
        enddo
    enddo

! fill cells at i = 1 and nz with i = 2 and nzm1

    do ni = nion1,nion2
        do k = 2,nlm1
            do j = 1,nf
                deni(1,j,k,ni)  = deni(2,j,k,ni)
                deni(nz,j,k,ni) = deni(nzm1,j,k,ni)
                ti(1,j,k,ni)    = ti(2,j,k,ni)
                ti(nz,j,k,ni)   = ti(nzm1,j,k,ni)
            enddo
        enddo
    enddo

    do k = 2,nlm1
        do j = 1,nf
            te(1,j,k)    = te(2,j,k)
            te(nz,j,k)   = te(nzm1,j,k)
        enddo
    enddo

!   fix at j = nf (interpolate)

!!$    do ni = nion1,nion2
!!$      do k = 2,nlm1
!!$      j = nf
!!$        do i = 1,nz
!!$          slope  = (blatp(i,j,k) - blatp(i,j-1,k)  ) /&
!!$                   (90.          - blatp(i,j-1,k))
!!$          deni(i,j,k,ni) = deni(i,j-1,k,ni) +  &
!!$                           slope * (deni_mnp(i,ni) - deni(i,j-1,k,ni)) 
!!$          deni(i,j,k,ni) = max(deni(i,j,k,ni),denmin)          
!!$        enddo
!!$      enddo
!!$    enddo

!   fix at j = nf (extapolate except at nz/2+1 - then interpolate)

!!$   j = nf
!!$
!!$    do ni = nion1,nion2
!!$      do k = 2,nlm1
!!$        do i = 1,nz/2
!!$          slope  = (blatp(i,j,k)   - blatp(i,j-1,k)  ) /&
!!$                   (blatp(i,j-1,k) - blatp(i,j-2,k))
!!$          deni(i,j,k,ni) = deni(i,j-1,k,ni) +  &
!!$                           slope * (deni(i,j-1,k,ni) - deni(i,j-2,k,ni)) 
!!$          deni(i,j,k,ni) = max(deni(i,j,k,ni),denmin)          
!!$        enddo
!!$      enddo
!!$    enddo
!!$
!!$    do ni = nion1,nion2
!!$      do k = 2,nlm1
!!$        do i = nz/2+2,nz
!!$          slope  = (blatp(i,j,k)   - blatp(i,j-1,k)  ) /&
!!$                   (blatp(i,j-1,k) - blatp(i,j-2,k))
!!$          deni(i,j,k,ni) = deni(i,j-1,k,ni) +  &
!!$                           slope * (deni(i,j-1,k,ni) - deni(i,j-2,k,ni)) 
!!$          deni(i,j,k,ni) = max(deni(i,j,k,ni),denmin)          
!!$        enddo
!!$      enddo
!!$    enddo
!!$
!!$    do ni = nion1,nion2
!!$      do k = 2,nlm1
!!$        i = nz/2+1
!!$          deni(i,j,k,ni) = 0.5 * ( deni(i-1,j,k,ni) + deni(i+1,j,k,ni) )
!!$          deni(i,j,k,ni) = max(deni(i,j,k,ni),denmin)          
!!$      enddo
!!$    enddo

    return
    end subroutine exb



!     ********************************************

!     VEXB_PHI

!     ********************************************

    subroutine vexb_phi (phi)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use exb_mod
    use grid_mod

    real :: phi(nnx,nny)
    real :: phihp(nfp1,nlp1)
    real :: phish(nfp1,nl)
    real :: phisp(nf,nlp1)

!      p face

    call phi_nfp1_nlp1(bradss,blonss,phihp,phi)

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                sinangx =  bdirpy(i,j,k) * ehpz(i,j,k) - &
                bdirpz(i,j,k) * ehpy(i,j,k)
                sinangy = -bdirpx(i,j,k) * ehpz(i,j,k) + &
                bdirpz(i,j,k) * ehpx(i,j,k)
                sinangz =  bdirpx(i,j,k) * ehpy(i,j,k) - &
                bdirpy(i,j,k) * ehpx(i,j,k)
                sinang  =  sqrt( sinangx * sinangx + &
                sinangy * sinangy + &
                sinangz * sinangz  )
                ehp(i,j,k) = -1. * &
                ( phihp(j,k+1) - phihp(j,k) )/ delhp(i,j,k) / &
                sinang
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                vexbp_phi(i,j,k) =  ehp(i,j,k) * sol / bmag / bmpf(i,j,k) &
                * tvexb0
            enddo
        enddo
    enddo

!      h face

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz
                sinangx =  bdirhy(i,j,k) * ephz(i,j,k) - &
                bdirhz(i,j,k) * ephy(i,j,k)
                sinangy = -bdirhx(i,j,k) * ephz(i,j,k) + &
                bdirhz(i,j,k) * ephx(i,j,k)
                sinangz =  bdirhx(i,j,k) * ephy(i,j,k) - &
                bdirhy(i,j,k) * ephx(i,j,k)
                sinang  =  sqrt( sinangx * sinangx + &
                sinangy * sinangy + &
                sinangz * sinangz  )
                eph(i,j,k) = -1. * &
                ( phihp(j+1,k) - phihp(j,k) )/ &
                delph(i,j,k)/sinang
            enddo
        enddo
    enddo

    do k = 1,nlp1
        do j = 1,nf-1
            do i = 1,nz
                vexbh_phi(i,j,k) =  -eph(i,j,k) * sol / bmag / bmhf(i,j,k) &
                * tvexb0
            enddo
        enddo
    enddo

!   fix at j = nf (extapolate except at nz/2+1 - then interpolate)

!!$   j = nf
!!$
!!$      do k = 1,nlp1
!!$        do i = 1,nz/2
!!$          slope  = (blatp(i,j,k)   - blatp(i,j-1,k)  ) /&
!!$                   (blatp(i,j-1,k) - blatp(i,j-2,k))
!!$          vexbh_phi(i,j,k) = vexbh_phi(i,j-1,k) +  &
!!$                           slope * (vexbh_phi(i,j-1,k) - vexbh_phi(i,j-2,k)) 
!!$        enddo
!!$      enddo
!!$
!!$      do k = 1,nlp1
!!$        do i = nz/2+2,nz
!!$          slope  = (blatp(i,j,k)   - blatp(i,j-1,k)  ) /&
!!$                   (blatp(i,j-1,k) - blatp(i,j-2,k))
!!$          vexbh_phi(i,j,k) = vexbh_phi(i,j-1,k) +  &
!!$                           slope * (vexbh_phi(i,j-1,k) - vexbh_phi(i,j-2,k)) 
!!$        enddo
!!$      enddo
!!$
!!$      do k = 1,nlp1
!!$        i = nz/2+1
!!$          vexbh_phi(i,j,k) = 0.5 * ( vexbh_phi(i-1,j,k) +  vexbh_phi(i+1,j,k) )
!!$      enddo

!      interpolate with vexbh_phi = 0 at pole

       j = nf
   
       do k = 1,nlp1
         do i = 1,nz 
           slope  = (blatp(i,j,k)   - blatp(i,j-1,k)  ) / &
                    (90.            - blatp(i,j-1,k))
           vexbh_phi(i,j,k) = vexbh_phi(i,j-1,k) * (1.- slope) &
             * tvexb0 
         enddo
       enddo

!      s face

    call phi_nfp1_nl  (bradsh,blonsh,phish,phi)
    call phi_nf_nlp1  (bradsp,blonsp,phisp,phi)

    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1
                eps(i,j,k) = -1. * &
                ( phish(j+1,k) - phish(j,k) )/delps(i,j,k)
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1
                ehs(i,j,k) = -1. * &
                ( phisp(j,k+1) - phisp(j,k) )/delhs(i,j,k)
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                vps = vexbp_phi(i,j,k) * &
                ( vpsnx(i,j,k) * xnorms(i,j,k) + &
                vpsny(i,j,k) * ynorms(i,j,k) + &
                vpsnz(i,j,k) * znorms(i,j,k)   )

                vph = vexbh_phi(i,j,k) * &
                ( vhsnx(i,j,k) * xnorms(i,j,k) + &
                vhsny(i,j,k) * ynorms(i,j,k) + &
                vhsnz(i,j,k) * znorms(i,j,k)   )
                            
                vexbs_phi(i,j,k) = (vps + vph) &
                * tvexb0
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nf
            vexbs_phi(nzp1,j,k) = vexbs_phi(nz,j,k)
        enddo
    enddo

    return
    end subroutine vexb_phi


!     ********************************************

!     current_jp_jphi

!     ********************************************

    subroutine current_jp_jphi

    use parameter_mod
    use variable_mod
    use conductance_mod
    use exb_mod
    use grid_mod

!     jp and jphi currents
!     gravity not included

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                ep        = eps(i,j,k)
                ephi      = ehs(i,j,k)
            !            ep        = eph(i,j,k)
            !            ephi      = ehp(i,j,k)
                boverc    = bmag * bms(i,j,k) / sol
                jp(i,j,k) = sigmap(i,j,k) * &
                ( ep + boverc * vnphi(i,j,k) ) &
                + sigmah(i,j,k) * &
                ( -ephi + boverc * vnp(i,j,k) )
                jphi(i,j,k) = sigmap(i,j,k) * &
                ( ephi - boverc * vnp(i,j,k) ) &
                + sigmah(i,j,k) * &
                ( ep + boverc * vnphi(i,j,k) )
            enddo
        enddo
    enddo


    return
    end subroutine current_jp_jphi


!     ********************************************

!     PHI_NFP1_NLP1

!     ********************************************


! takes a variable in one grid (xu,zu):(nf,nl)
! and interpolates to another grid (x0,z0):(nx,ny)
! the nf,nf grid need not be orthogonal

    subroutine phi_nfp1_nlp1(brad,blon,phiout,phiin)

    use parameter_mod
    use variable_mod
    use message_passing_mod
    use grid_mod

    real :: brad(nfp1,nlp1),blon(nfp1,nlp1),phiout(nfp1,nlp1)
    real :: phiin(nnx,nny)

    do j = 1,nlp1
        do i = 1,nfp1
            phiout(i,j) = 0.
        enddo
    enddo

    if ( taskid /= 1 .AND. taskid /= numworkers ) then
        do k = 1,nlp1
            do j = nf,2,-1
                kk = (taskid-1)*(nl-2) + (k-1)
                jj = j -1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == 1 ) then
        do k = 1,nlp1
            do j = nf,2,-1
                kk = k - 1
                if ( k == 1 ) kk = nnx - 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == numworkers ) then
        do k = 1,nlp1
            do j = nf,2,-1
                kk = (taskid - 1)*(nl-2) + (k-1)
                if ( k == nlp1   ) kk = 2
                if ( k == nlp1-1 ) kk = 1
                jj = j -1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    do k = 1,nlp1
        phiout(1,k)     = phiout(2,k)
        slope           = (blatp(nz-1,nf+1,k) - blatp(nz-1,nf-1,k))/ &
        (blatp(nz-1,nf,k)   - blatp(nz-1,nf-1,k))
        phiout(nf+1,k)  = phiout(nf-1,k) + &
        (phiout(nf,k) - phiout(nf-1,k)) * slope
         

    enddo

    return
    end subroutine phi_nfp1_nlp1

!     ********************************************

!     PHI_NFP1_NL

!     ********************************************


! takes a variable in one grid (xu,zu):(nf,nl)
! and interpolates to another grid (x0,z0):(nx,ny)
! the nf,nf grid need not be orthogonal

    subroutine phi_nfp1_nl(brad,blon,phiout,phiin)

    use parameter_mod
    use variable_mod
    use message_passing_mod
    use grid_mod

    real :: brad(nfp1,nl),blon(nfp1,nl)
    real :: phiout(nfp1,nl)
    real :: phiin(nnx,nny)

    do j = 1,nl
        do i = 1,nfp1
            phiout(i,j) = 0.
        enddo
    enddo

    if ( taskid /= 1 .AND. taskid /= numworkers ) then
        do k = 1,nl
            do j = nf,2,-1
                kk = (taskid-1)*(nl-2) + (k-1)
                jj = j -1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == 1 ) then
        do k = 1,nl
            do j = nf,2,-1
                kk = k - 1
                if ( k == 1 ) kk = nnx - 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == numworkers ) then
        do k = 1,nl
            do j = nf,2,-1
                kk = (taskid - 1)*(nl-2) + (k-1)
                if ( k == nl   ) kk = 1
                jj = j -1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    do k = 1,nl
        phiout(1,k)     = phiout(2,k)
        slope           = (blatp(nz-1,nf+1,k) - blatp(nz-1,nf-1,k))/ &
        (blatp(nz-1,nf,k)   - blatp(nz-1,nf-1,k))
        phiout(nf+1,k)  = phiout(nf-1,k) + &
        (phiout(nf,k) - phiout(nf-1,k)) * slope
         

    enddo

    return
    end subroutine phi_nfp1_nl



!     ********************************************

!     PHI_NF_NLP1

!     ********************************************


    subroutine phi_nf_nlp1(brad,blon,phiout,phiin)

    use parameter_mod
    use variable_mod
    use message_passing_mod
    use grid_mod

    real :: brad(nf,nlp1),blon(nf,nlp1),phiout(nf,nlp1)
    real :: phiin(nnx,nny)

    do j = 1,nlp1
        do i = 1,nf
            phiout(i,j) = 0.
        enddo
    enddo

    if ( taskid /= 1 .AND. taskid /= numworkers ) then
        do k = 1,nlp1
            do j = nf-1,2,-1
                kk = (taskid-1)*(nl-2) + (k-1)
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == 1 ) then
        do k = 1,nlp1
            do j = nf-1,2,-1
                kk = k - 1
                if ( k == 1 ) kk = nnx - 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == numworkers ) then
        do k = 1,nlp1
            do j = nf-1,2,-1
                kk = (taskid - 1)*(nl-2) + (k-1)
                if ( k == nlp1   ) kk = 2
                if ( k == nlp1-1 ) kk = 1
                jj = j -1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    do k = 1,nlp1
        phiout(nf,k)  = phiout(nf-1,k) * &
        (blatp(nz-1,nf,k)   - 90.) / &
        (blatp(nz-1,nf-1,k) - 90.)
        phiout(1,k)   = phiout(2,k)
    enddo

    return
    end subroutine phi_nf_nlp1

