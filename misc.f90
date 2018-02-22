
!******************************************
!******************************************

!             courant

!******************************************
!******************************************

    subroutine courant

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use exb_mod
    use grid_mod

! parallel motion

    dtnew = 1.e6
    do k = nion1,nion2
        do l = 1,nl
            do j = 1,nf
                do i = 1,nz
                    dt1 = dels(i,j,l) / amax1(1.,abs(vsi(i,j,l,k)))
                    if ( dt1 <= dtnew ) then
                        dtnew = dt1
                        i0    = i
                        j0    = j
                        k0    = k
                        l0    = l
                    endif
                enddo
            enddo
        enddo
    enddo

! perpendicular motion

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                dts = xdels(i,j,k) / amax1(1.,abs(vexbs(i,j,k)))
                dtp = xdelp(i,j,k) / amax1(1.,abs(vexbp(i,j,k)))
                dth = xdelh(i,j,k) / amax1(1.,abs(vexbh(i,j,k)))
                dt1 = amin1 ( dts,dtp,dth )
                if ( dt1 <= dtnew ) then
                    dtnew = dt1
                    i0    = i
                    j0    = j
                    k0    = k
                endif
            enddo
        enddo
    enddo

    if ( dtnew <= .01 ) then
        print *,' Time step too small'
        stop
    elseif ( dtnew >= 5e4 ) then
        print *,' Time step too big: dtnew',dtnew
        stop
    endif
    dtnew = 1.5 * dtnew
    if ( dtnew <= dt      ) dt = dtnew
    if ( dtnew > dt*1.2  ) dt = dt * 1.2

    return
    end subroutine courant


! *********************

!     smoothz

! *********************

    subroutine smoothz(finout,ncomp)
     
    use parameter_mod
     
    dimension finout(nz), tempz(nz)
     

! This is the binomial filter (in x space) as described in
! Birdsall appendix C.
! We have the choice of a compensating filter or not.
! if ncomp=0, no compensation, else compensation

     
! do smoothz in the z direction
     
    do i = 1,nz
        ip1 = i +1
        if(i == nz) ip1 = 1
        im1 = i -1
        if(i == 1) im1 = nz
        tempz(i) = .25*(finout(im1) +2.*finout(i) &
        +finout(ip1))
    enddo
    do i = 1,nz
        finout(i) = tempz(i)
    enddo
     
    if ( ncomp /= 0 ) then
         
    ! put in compensator
    ! the alogrithm below is equivalent to
    ! fftmp(i)=(1./16.)*(-ff0(i-2)+4.*ff0(i-1)+10.*ff0(i)+4.*ff0(i+1)-ff0(i+2))
         
    ! do compensation in the z direction
         
        const = sqrt(1.4571072)
        do i = 1,nz
            ip1 = i +1
            if(i == nz) ip1 = 1
            finout(i) = const*(finout(i) -.171573*finout(ip1))
        enddo
        do i = nz,1,-1
            im1 = i -1
            if(i == 1) im1 = nz
            finout(i) = const*(finout(i) -.171573*finout(im1))
        enddo
         
    endif
     
    return
    end subroutine smoothz



!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************

    subroutine smoothx(f)
     
    use parameter_mod

    parameter ( nnxp2 = nnx + 2, nnyp2 = nny + 2 )
    parameter ( nnxp1 = nnx + 1, nnyp1 = nny + 1 )
             
    real :: f(nnx,nny),f0(nnx+2,nny+2)
     
    u12 = 1.
     
    do j = 1,nny
        do i = 1,nnx
            f0(i+1,j+1) = f(i,j)
        enddo
    enddo
     
!       zero-gradient in x
     
!        do j = 2,nnyp1
!          f0(1,j)     = f0(2,j)
!          f0(nnx+2,j)  = f0(nnx+1,j)
!        enddo

!       periodic in x
     
    do j = 2,nnyp1
        f0(1,j)     = f0(nnx,j)
        f0(nnx+2,j)  = f0(3,j)
    enddo
     
!       zero gradient in y
     
    do i = 1,nnxp2
        f0(i,1)    = f0(i,2)
        f0(i,nnyp2) = f0(i,nnyp1)
    enddo
     
!       sweep in x (1/2)
        
    do j = 2,nnyp1
        do i = 1,nnxp1
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i+1,j) )
        enddo
    enddo
     
    do j = 2,nnyp1
        do i = nnxp2,2,-1
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i-1,j) )
        enddo
    enddo
     
!       now get f
     
    do j = 1,nny
        do i = 1,nnx
            f(i,j)  = f0(i+1,j+1)
        enddo
    enddo
     
    return
    end subroutine smoothx

!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************


    subroutine smoothy(f)
     
    use parameter_mod

    parameter ( nnxp2 = nnx + 2, nnyp2 = nny + 2 )
    parameter ( nnxp1 = nnx + 1, nnyp1 = nny + 1 )
             
    real :: f(nnx,nny),f0(nnx+2,nny+2)
     
    u12 = 1.
     
    do j = 1,nny
        do i = 1,nnx
            f0(i+1,j+1) = f(i,j)
        enddo
    enddo
            
!       zero-gradient in x
     
!        do j = 2,nnyp1
!          f0(1,j)     = f0(2,j)
!          f0(nnx+2,j)  = f0(nnx+1,j)
!        enddo
     
!       periodic in x
     
    do j = 2,nnyp1
        f0(1,j)     = f0(nnx,j)
        f0(nnx+2,j)  = f0(3,j)
    enddo

!       zero gradient in y
     
    do i = 1,nnxp2
        f0(i,1)    = f0(i,2)
        f0(i,nnyp2) = f0(i,nnyp1)
    enddo
     
!       sweep in y (1/2)
     
    do j = 1,nnyp1
        do i = 1,nnxp2
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i,j+1) )
        enddo
    enddo
     
    do j = nnyp2,2,-1
        do i = 1,nnxp2
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i,j-1) )
        enddo
    enddo
     
!       now get f
     
    do j = 1,nny
        do i = 1,nnx
            f(i,j)  = f0(i+1,j+1)
        enddo
    enddo
     
    return
    end subroutine smoothy
     



!******************************************
!******************************************
!          splinenr
!    (from numerical recipes)
!******************************************
!******************************************

    subroutine splinenr(x,y,n,yp1,ypn,y2)
    parameter (nmax=200)
    dimension x(n),y(n),y2(n),u(nmax)
    if (yp1 > .99e30) then
        y2(1)=0.
        u(1)=0.
    else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    11 END DO
    if (ypn > .99e30) then
        qn=0.
        un=0.
    else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    12 END DO
    return
    end subroutine splinenr

!******************************************
!******************************************
!          splintnr
!    (from numerical recipes)
!******************************************
!******************************************

    subroutine splintnr(xa,ya,y2a,n,x,y)
    dimension xa(n),ya(n),y2a(n)
    klo=1
    khi=n
    1 if (khi-klo > 1) then
        k=(khi+klo)/2
        if(xa(k) > x)then
            khi=k
        else
            klo=k
        endif
        goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h == 0.) print *, 'bad xa input. from splintnr'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+ &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    return
    end subroutine splintnr


