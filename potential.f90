!       **********************************
!       **********************************

!             POTPPHI

!       **********************************
!       **********************************

    subroutine potpphi(phi,dphi,hrut)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use conductance_mod
    use misc_mod
    use grid_mod

!        parameter ( nyextra = 10, nnyt = nny + nyextra )

    real :: hipcp_pot(nnx,nnyt),hipcphi_pot(nnx,nnyt)
    real :: hidphig_pot(nnx,nnyt),hidphiv_pot(nnx,nnyt)
    real :: hidpg_pot(nnx,nnyt),hidpv_pot(nnx,nnyt)
    real :: hihcm_pot(nnx,nnyt)
    real :: hipc_pot(nnx,nnyt)
    real :: hihc_pot(nnx,nnyt)
    real :: hidv_pot(nnx,nnyt)

    real :: phi(nnx,nny)
    real :: dpreal(nnyt),preal(nnyt)
    real :: dbang(nnx,nnyt),blang(nnx)
    real*8 :: dphireal(nnx+1)
    real*8 :: dxij,dxip1j
    real*8 :: f11_lb(nnx+1),f11_ub(nnx+1)

    real*8 :: dphi(nnx+1,nnyt),si(nnx+1,nnyt),sih(nnx+1,nny)
    real*8 :: a1(nnx+1,nnyt),a2(nnx+1,nnyt),a3(nnx+1,nnyt)
    real*8 :: a4(nnx+1,nnyt),a5(nnx+1,nnyt)
    real*8 :: dphi0(nnx+1,nnyt)

    real :: ylonp(nnx),ylatp(nny)
    real :: zigm11(nnx,nny),zigm22(nnx,nny),zigm2(nnx,nny)
    real :: rim1(nnx,nny),rim2(nnx,nny)

!    real :: vexbpphi(nnx,nny),vexbhphi(nnx,nny)
    real :: phivs(nnx,nnyt)

    real :: phi_weimer(nnx,nny)

    data idpinter / 1 /
    data ipcrit   / 1 /

    if ( idpinter == 1 ) then

        nzh  = nz / 2
              
        do j = 1,nny
            jj         = j + nyextra
            dpreal(jj) = ppt(nzh,j+1,1) - ppt(nzh,j,1)
            preal(jj)  = ppt(nzh,j,1)
!            print *,'in pot',jj,preal(jj)
        enddo

        do j = nyextra,1,-1
            dpreal(j) = dpreal(j+1) * 1.4
            preal(j)  = preal(j+1) - dpreal(j)
        enddo

        do i = 1,nnx+1
            dphireal(i) = (blonp0t(i+1)-blonp0t(i))*pie/180.
        enddo

    ! define blang
    ! dbang (new dy)

        do i = 1,nnx-1
            do j = 1,nny
                jj          = j + nyextra
                blang(i)    = -blonpt(nz/2+1,j,i) * pie / 180.
            enddo
        enddo

        blang(nnx) = blang(1)

    !        initialize dphi and dphi0

        do j = 1,nnyt
            do i = 1,nnx+1
                dphi (i,j) = 0.
                dphi0(i,j) = 0.
            enddo
        enddo

        if (restart) then
            open(1232,file='dphi.rst',form='unformatted')
            read(1232) dphi0
            close(1232)
        endif

        idpinter = 0

    endif

!       set up conductances and driver for potential equation
!       zero-gradient in phi (x); zero-gradient in p (y)
!       note: transpose variables

    do j = 1,nny
        jj   = j + nyextra
        do i = 2,nnx-1
            hipcp_pot(i,jj)    = 0.25 * ( hipcpt(j,i-1)   + &
            hipcpt(j,i)     + &
            hipcpt(j+1,i-1) + &
            hipcpt(j+1,i)     )
            hihcm_pot(i,jj)    = 0.25 * ( hihcmt(j,i-1)   + &
            hihcmt(j,i)     + &
            hihcmt(j+1,i-1) + &
            hihcmt(j+1,i)     )
            hipcphi_pot(i,jj)  = 0.25 * ( hipcphit(j,i-1)   + &
            hipcphit(j,i)     + &
            hipcphit(j+1,i-1) + &
            hipcphit(j+1,i)     )
            hidphig_pot(i,jj)  = 0.25 * ( hidphigt(j,i-1)   + &
            hidphigt(j,i)     + &
            hidphigt(j+1,i-1) + &
            hidphigt(j+1,i)     )
            hidpg_pot(i,jj)    = 0.25 * ( hidpgt(j,i-1)   + &
            hidpgt(j,i)     + &
            hidpgt(j+1,i-1) + &
            hidpgt(j+1,i)     )
            hidphiv_pot(i,jj)  = 0.25 * ( hidphivt(j,i-1)   + &
            hidphivt(j,i)     + &
            hidphivt(j+1,i-1) + &
            hidphivt(j+1,i)     )
            hidpv_pot(i,jj)    = 0.25 * ( hidpvt(j,i-1)   + &
            hidpvt(j,i)     + &
            hidpvt(j+1,i-1) + &
            hidpvt(j+1,i)     )
            hipc_pot(i,jj)     = 0.25 * ( hipct(j,i-1)   + &
            hipct(j,i)     + &
            hipct(j+1,i-1) + &
            hipct(j+1,i)     )
            hihc_pot(i,jj)     = 0.25 * ( hihct(j,i-1)   + &
            hihct(j,i)     + &
            hihct(j+1,i-1) + &
            hihct(j+1,i)     )
            hidv_pot(i,jj)     = 0.25 * ( hidvt(j,i-1)   + &
            hidvt(j,i)     + &
            hidvt(j+1,i-1) + &
            hidvt(j+1,i)     )
        enddo
    enddo

    do j = nyextra+1,nnyt
        hipcp_pot(1,j)   = 0.5 * (hipcp_pot(2,j)    + &
        hipcp_pot(nnx-1,j) )
        hihcm_pot(1,j)   = 0.5 * (hihcm_pot(2,j)    + &
        hihcm_pot(nnx-1,j) )
        hipcphi_pot(1,j) = 0.5 * (hipcphi_pot(2,j)  + &
        hipcphi_pot(nnx-1,j) )
        hidphig_pot(1,j) = 0.5 * ( hidphig_pot(2,j) + &
        hidphig_pot(nnx-1,j) )
        hidpg_pot(1,j)   = 0.5 * ( hidpg_pot(2,j) + &
        hidpg_pot(nnx-1,j) )
        hidphiv_pot(1,j) = 0.5 * ( hidphiv_pot(2,j) + &
        hidphiv_pot(nnx-1,j) )
        hidpv_pot(1,j)   = 0.5 * ( hidpv_pot(2,j) + &
        hidpv_pot(nnx-1,j) )
        hipc_pot(1,j)    = 0.5 * (hipc_pot(2,j)    + &
        hipc_pot(nnx-1,j) )
        hihc_pot(1,j)    = 0.5 * (hihc_pot(2,j)    + &
        hihc_pot(nnx-1,j) )
        hidv_pot(1,j)    = 0.5 * (hidv_pot(2,j)    + &
        hidv_pot(nnx-1,j) )
    enddo

    do j = nyextra+1,nnyt
        hipcp_pot(nnx,j)   = hipcp_pot(1,j)
        hihcm_pot(nnx,j)   = hihcm_pot(1,j)
        hipcphi_pot(nnx,j) = hipcphi_pot(1,j)
        hidphig_pot(nnx,j) = hidphig_pot(1,j)
        hidpg_pot(nnx,j)   = hidpg_pot(1,j)
        hidphiv_pot(nnx,j) = hidphiv_pot(1,j)
        hidpv_pot(nnx,j)   = hidpv_pot(1,j)
        hipc_pot(nnx,j)    = hipc_pot(1,j)
        hihc_pot(nnx,j)    = hihc_pot(1,j)
        hidv_pot(nnx,j)    = hidv_pot(1,j)
    enddo

    do j = nyextra,1,-1
        do i = 1,nnx
            hipcp_pot(i,j)   = hipcp_pot(i,j+1)   * .02
            hihcm_pot(i,j)   = hihcm_pot(i,j+1)   * .02
            hipcphi_pot(i,j) = hipcphi_pot(i,j+1) * .02
            hidphig_pot(i,j) = hidphig_pot(i,j+1) * .01
            hidpg_pot(i,j)   = hidpg_pot(i,j+1)   * .01
            hidphiv_pot(i,j) = hidphiv_pot(i,j+1) * .01
            hidpv_pot(i,j)   = hidpv_pot(i,j+1)   * .01
            hipc_pot(i,j)    = hipc_pot(i,j+1)    * .02
            hihc_pot(i,j)    = hihc_pot(i,j+1)    * .02
            hidv_pot(i,j)    = hidv_pot(i,j+1)    * .01
        enddo
    enddo

    k   = 0
    k0  = 10

     i   = nnx/2
     im1 = nnx/2-1
     j   = nnyt/1

    do while ( k <= k0 )

    !     div j = 0 differencing (for Pedersen and Hall via iteration)

        do j = 2,nnyt-1
            do i = 1,nnx

                im1  = i - 1
                ip1  = i + 1
                jm1  = j - 1
                jp1  = j + 1

                dxij   = dphireal(i)
                dxip1j = dphireal(ip1)

                dyij    = dpreal(j)
                dyijp1  = dpreal(jp1)

                if ( i == 1   ) im1 = nnx - 1
                if ( i == nnx ) ip1 = 2

                delx    = dxij + dxip1j
                dely    = dyij + dyijp1
                 
                delx_inv  = 1. / delx
                dely_inv  = 1. / dely
                delxy_inv = delx_inv * dely_inv

!!$                pcphimx  = (0.5 * ( hipcphi_pot(im1,j) &
!!$                + hipcphi_pot(i,j)   ))
!!$                pcphipx  = (0.5 * ( hipcphi_pot(i,j) &
!!$                + hipcphi_pot(ip1,j) ))
!!$
!!$                pcpmy  = (0.5 * ( hipcp_pot(i,jm1) + hipcp_pot(i,j)   ))
!!$                pcppy  = (0.5 * ( hipcp_pot(i,j)   + hipcp_pot(i,jp1) ))

          pcphimx  = (0.5 * ( hipcphi_pot(im1,j) &
                     + hipcphi_pot(i,j)   )) / preal(j)
          pcphipx  = (0.5 * ( hipcphi_pot(i,j)   &
                     + hipcphi_pot(ip1,j) )) / preal(j)

! below are 2 ways of writing pcpmy/pcppy
! both give basically the same vertical drifts and phi
! (at least after 10 time steps)

!!$          pcpmy  = (0.5 * ( hipcp_pot(i,jm1) + hipcp_pot(i,j)   )) &
!!$                 * (0.5 * ( preal(jm1) + preal(j) ) )  
!!$          pcppy  = (0.5 * ( hipcp_pot(i,j)   + hipcp_pot(i,jp1) )) &
!!$                 * (0.5 * ( preal(j) + preal(jp1) ) )

          pcpmy  = 0.5 * ( hipcp_pot(i,jm1) * preal(jm1) + &
                           hipcp_pot(i,j)   * preal(j)   ) 
          pcppy  = 0.5 * ( hipcp_pot(i,j)   * preal(j)   + &
                           hipcp_pot(i,jp1) * preal(jp1) ) 

                if (hall) then

                    hcmjp = 0.5 * ( hihcm_pot(i,jp1) + hihcm_pot(i,j)   )
                    hcmjm = 0.5 * ( hihcm_pot(i,j)   + hihcm_pot(i,jm1) )
                    dhcy  = hcmjp - hcmjm

                    hcmip = 0.5 * ( hihcm_pot(i,j)   + hihcm_pot(ip1,j) )
                    hcmim = 0.5 * ( hihcm_pot(im1,j) + hihcm_pot(i,j)   )
                    dhcx  = hcmip - hcmim

                else

                    dhcy   = 0.
                    dhcx   = 0.
                    hcmjp  = 0.
                    hcmjm  = 0.
                    hcmip  = 0.
                    hcmim  = 0.

                endif

                fphivmx   = 0.5 * ( hidphiv_pot(im1,j) + hidphiv_pot(i,j)   )
                fphivpx   = 0.5 * ( hidphiv_pot(i,j)   + hidphiv_pot(ip1,j) )
                dfphivx   = fphivpx - fphivmx

                fphigmx   = 0.5 * ( hidphig_pot(im1,j) + hidphig_pot(i,j)   )
                fphigpx   = 0.5 * ( hidphig_pot(i,j)   + hidphig_pot(ip1,j) )
                dfphigx   = fphigpx - fphigmx


                fpvmy   = 0.5 * ( hidpv_pot(i,jm1) + hidpv_pot(i,j)   )
                fpvpy   = 0.5 * ( hidpv_pot(i,j)   + hidpv_pot(i,jp1) )
                dfpvy   = fpvpy - fpvmy

                fpgmy   = 0.5 * ( hidpg_pot(i,jm1) + hidpg_pot(i,j)   )
                fpgpy   = 0.5 * ( hidpg_pot(i,j)   + hidpg_pot(i,jp1) )
                dfpgy   = fpgpy - fpgmy

                a11    = 2. * delx_inv * pcphimx / dxij
                a12    = delxy_inv * dhcy

                a1(i,j) = a11 + a12

                a21     = 2. * dely_inv * pcpmy / dyij
                a22     = delxy_inv * dhcx

                a2(i,j) = a21 - a22

                a41     = 2. * dely_inv * pcppy / dyijp1
                a42     = delxy_inv * dhcx

                a4(i,j) = a41 + a42

                a51     = 2. * delx_inv * pcphipx / dxip1j
                a52     = delxy_inv * dhcy

                a5(i,j) = a51 - a52

                a3(i,j) = -a51 - a11 - a41 - a21

                sip1jp1 = -delxy_inv * (hcmip - hcmjp) * dphi0(ip1,jp1)

                sip1jm1 =  delxy_inv * (hcmip - hcmjm) * dphi0(ip1,jm1)

                sim1jp1 =  delxy_inv * (hcmim - hcmjp) * dphi0(im1,jp1)

                sim1jm1 = -delxy_inv * (hcmim - hcmjm) * dphi0(im1,jm1)

                si(i,j) = &
                  2. * dfphigx * delx_inv &
                + 2. * dfphivx * delx_inv &
                - 2. * dfpgy   * dely_inv &
                + 2. * dfpvy   * dely_inv &
                + sip1jp1 + sip1jm1 + sim1jp1 + sim1jm1

!!$        if ( i == 48 .and. j == 20 ) then 
!!$          print *,'s1',dfphigx,dfphivx,delx_inv
!!$          print *,'s2',dfpgy,dfpvy,dely_inv
!!$        endif

            enddo
        enddo

    !     zero gradient in y

        do i = 1,nnx
            a1(i,1)   = a1(i,2)
            a2(i,1)   = a2(i,2)
            a3(i,1)   = a3(i,2)
            a4(i,1)   = a4(i,2)
            a5(i,1)   = a5(i,2)
            si(i,1)   = si(i,2)

            dpp          = dpreal(nnyt) / dpreal(nnyt-1)

            a1(i,nnyt)   = dpp*(a1(i,nnyt-1)-a1(i,nnyt-2))+a1(i,nnyt-1)
            a2(i,nnyt)   = dpp*(a2(i,nnyt-1)-a2(i,nnyt-2))+a2(i,nnyt-1)
            a3(i,nnyt)   = dpp*(a3(i,nnyt-1)-a3(i,nnyt-2))+a3(i,nnyt-1)
            a4(i,nnyt)   = dpp*(a4(i,nnyt-1)-a4(i,nnyt-2))+a4(i,nnyt-1)
            a5(i,nnyt)   = dpp*(a5(i,nnyt-1)-a5(i,nnyt-2))+a5(i,nnyt-1)
            si(i,nnyt)   = dpp*(si(i,nnyt-1)-si(i,nnyt-2))+si(i,nnyt-1)

        enddo

        do j = 1,nnyt
            a1(nnx+1,j) = a1(2,j)
            a2(nnx+1,j) = a2(2,j)
            a3(nnx+1,j) = a3(2,j)
            a4(nnx+1,j) = a4(2,j)
            a5(nnx+1,j) = a5(2,j)
            si(nnx+1,j) = si(2,j)
        enddo

        do i = 2,nnx
            f11_lb(i)  = 0.
        enddo

        f11_lb(1)     = f11_lb(nnx)
        f11_lb(nnx+1) = f11_lb(2)

        phi90 = 0.

        do i = 2,nnx
!               f11_ub(i)  = 0.

            dbangi90     = (blatpt(nz-1,nf-2,1)-90.)
            dbangif      = (blatpt(nz-1,nf-2,1) - blatpt(nz-1,nf-1,1))
            dbangf90     = (blatpt(nz-1,nf-1,1)-90.)

            dphi(i,nnyt) = ( dphi(i,nnyt-1) * dbangf90 + &
            phi90          * dbangif  ) / dbangi90

        enddo


        f11_ub(1)     = f11_ub(nnx)
        f11_ub(nnx+1) = f11_ub(2)

        if ( lmadala ) then
            call madala_sevp(a1,a2,a5,a4,a3,si,dphi,f11_lb,f11_ub)
        else
            do j = 1,nnyt
                do i = 1,nnx+1
                    dphi(i,j) = 0.
                enddo
            enddo
        endif

        i00 = nnx/2
        j00 = nnyt/2

        err0 = abs(dphi0(i00,j00)-dphi(i00,j00))/ &
        abs(dphi0(i00,j00)+1.e-5)
        k = k + 1
        print *,k,err0
         
        tol_phi  = 2.e-2

        if ( err0 <= tol_phi ) k0 = -10

        do j = 1,nnyt
            do i = 1,nnx+1
                dphi0(i,j) = dphi(i,j)
            enddo
        enddo

    enddo

!     parameters for volland/stern/macilwane potential

    fkp  = 6.

    akp_coef_i = 11.25e-2
    akp_coef_f = 45.00

    if ( hrut < storm_ti ) akp_coef = akp_coef_i
    if ( hrut > storm_tf ) akp_coef = akp_coef_f
    if ( hrut >= storm_ti .AND. hrut <= storm_tf ) &
    akp_coef = akp_coef_i + &
    (akp_coef_f - akp_coef_i) * &
    (hrut - storm_ti)/(storm_tf - storm_ti)

    akp  = akp_coef / ( (1.-0.159*fkp+0.0093*fkp*fkp) ** 3 )

!     only for non-corotating code

!     rotate phi so that potential is
!     aligned with local time midnight/noon
!       angut  in degrees
!       angutr in radians

    angrot = 360. - plon * 180. / pie
    hr24   = mod (hrut,24.)
    angut  = hr24 * 2. * pie / 24. * rtod - angrot
    if ( angut > 360. ) angut = angut - 360.
    if ( angut < 0.   ) angut = angut + 360.
    angutr = angut * pie / 180.

    if ( lcr ) then
        angutr  = 0.
        angut   = 0.
    endif

    if ( lweimer ) then
        call weimer(phi_weimer,angut,hrut)
    else
        do j = nyextra+1,nnyt
            do i = 1,nnx
                jj = j - nyextra
                phi_weimer(i,jj) = 0.
            enddo
        enddo
    endif

    do j = nyextra+1,nnyt

        do i = 1,nnx
            jj = j - nyextra
            phivs(i,jj)   = 0.
            if ( lvs ) &
            phivs(i,jj)   = -akp * preal(j) * preal(j) * &
            sin(blang(i)-angutr) &
            / 300.
            phicorot    = 0.
            if ( lcr ) &
            phicorot    = -92.e3 / 300. / preal(j) * bmag /.31
            phi(i,jj)   = dphi(i,j) + phicorot + phivs(i,jj) + &
            phi_weimer(i,jj)
        enddo
    enddo

    print *,phi(40,35)


    return

    end subroutine potpphi


!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************

    subroutine weimer(phi_weimer,angut,hrut)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use misc_mod
    use grid_mod

    real :: phi_weimer_real(nfp1,nlt+1)
    real :: phi_weimer_interp(nfp1,nlt+1)
    real :: phi_weimer(nnx,nny)
    real :: hrutw1,hrutw2

    data iread_weimer / 1 /

!       read in weimer at first time

    if ( iread_weimer == 1 ) then
        open(810,file='phi_weimer.inp',form='unformatted')
        read(810) hrutw1
        if ( restart ) then
            open (881,file='nweimer.rst',form='formatted')
            read(881,*) nweimer
            print *,'nweimer',nweimer
            do i = 1,nweimer-1
                read(810) phi_weimer_real
                read(810) hrutw2
            enddo
        else
            read(810) phi_weimer_real
            read(810) hrutw2
            print *,'hrutw2 = ',hrut,hrutw2
            nweimer       = 1
            print *,'nweimer',nweimer
        endif
        iread_weimer  = 0
    endif

!       read in weimer for subsequent time steps

    if ( hrut >= hrutw2 ) then
        read(810) phi_weimer_real
        read(810) hrutw2
        print *,'hrutw2 = ',hrut,hrutw2
        nweimer         = nweimer + 1
    endif

!       save nweimer for restart

    open (881,file='nweimer.rst',form='formatted')
    write(881,*) nweimer
    close(881)

!       requires sami3 and weimer potential to
!       have the same latitude
!       (from sami3_rcm interp_rcm subroutine)

    do j = 1,nny
        do i = 1,nnx
            phi_weimer(i,j) = 0.
        enddo
    enddo

    dlon = 360. / float(nlt)

    do k = 1,nlt+1
        do j = nfp1,1,-1
            thlon = mod(blonp0t(k+1) + angut,360.)
            if ( thlon < 0. ) thlon = thlon + 360.
            nlon = thlon/dlon + 1
            dnlon = thlon/dlon - nlon + 1
            phi_weimer_interp(j,k) = &
            + dnlon * phi_weimer_real(j,nlon+1) &
            + (1. - dnlon) * phi_weimer_real(j,nlon)
        !            if (j.eq.121) then
        !              print *,'k121',k,thlon,nlon,dnlon,
        !     .                  phi_weimer_real(j,nlon),phi_weimer_interp(j,k)
        !            endif
        !            if (j.eq.122) then
        !              print *,'k122',k,thlon,nlon,dnlon,
        !     .                  phi_weimer_real(j,nlon),phi_weimer_interp(j,k)
        !            endif
        enddo
    enddo

!       here nnx = nlt + 1
!            nny = nf  - 1

    do j = 1,nny
        do i = 1,nnx
            phi_weimer(i,j) = phi_weimer_interp(j,i)
        enddo
    enddo

!       try to fix weimer potential problem with
!       simple linear interpolation

    j    = nny - 1
    do i = 1,nnx
        phi_weimer(i,j) = 0.5 * ( phi_weimer(i,j-1) + &
        phi_weimer(i,j+1)   )
    enddo

    return
    end subroutine weimer
