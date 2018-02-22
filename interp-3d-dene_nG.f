
!
!      INTERP-3D-DENE.F
!

!      takes deneu in sami3 coordinates (zaltu,glatu,glonu):(nz,nf,nl,nt)
!      and interpolates to a regular grid (zalt0,glat0,glon0):(nx,ny,nl,nt)
!        - the grid in the z-direction can be nonuniform (gamy)

       include "param_diag.inc"

!      arrays of 'real' data

       real glatu(nz,nf,nl),zaltu(nz,nf,nl),glonu(nz,nf,nl)
       real deneu(nz,nf,nl,nt),denetmp(nz,nf,nl)
       real xu(nz,nf,nl),zu(nz,nf,nl)
       real zualtmax(nf,nl)    !gj


!      arrays of interpolated data

       real    x0(nx,ny,nl),z0(nx,ny,nl)
       integer i0(nx,ny,nl),j0(nx,ny,nl)

       real glat0(nx,ny,nl),glon0(nx,ny,nl),zalt0(nx,ny,nl)
       real dene0(nx,ny,nl,nt)

!      open 'real' data files

       open(unit= 9,file='glonu.dat', form='unformatted')
       open(unit=12,file='deneu.dat' , form='unformatted')

!      open map file

       open(unit=60,file='map.dat', form='unformatted')

!      open interpolated data files

       open(unit=31,file='glat0.dat', form='unformatted')
       open(unit=32,file='glon0.dat', form='unformatted')
       open(unit=33,file='zalt0.dat', form='unformatted')
       open(unit=34,file='dene0.dat', form='unformatted')

       open(unit=45,file='x0.dat', form='unformatted')
       open(unit=46,file='z0.dat', form='unformatted')
       open(unit=55,file='xu.dat', form='unformatted')
       open(unit=56,file='zu.dat', form='unformatted')

!      added smoother to dene0: ipsmooth  
!      ipsmooth = 1 or 2 seems best 

!      smoothing parameter 
 
       ipsmooth  = 2 
 

!      read in 'real' data

       read( 9) glonu

       read(45) x0
       read(46) z0
       read(55) xu
       read(56) zu

       do l = 1,nt
         print *,' reading data: l = ',l
         read(12) denetmp
         do k = 1,nl
           do j = 1,nf
             do i = 1,nz
               deneu(i,j,k,l) = denetmp(i,j,k)
             enddo
           enddo
         enddo
       enddo

       print *,'getting interior points from map.dat'

       read(60) i0,j0

!      set default value of interpolated data
!      (in this case the electron density)

       print *,'initializing dene0'

       do n = 1,nt
         do k = 1,nl
           do j = 1,ny
             do i = 1,nx
                dene0(i,j,k,n) = 1.
             enddo
           enddo
         enddo
       enddo

!      lay down interpolated data on interpolated grid
!      using an area weighted lay down scheme

       print *,'setting dene0'

       do n = 1,nt
         do k = 1,nl
           do iy = 1,ny
             do ix = 1,nx
                i           = i0(ix,iy,k)
                j           = j0(ix,iy,k)
                if ( i .ne. 0 .and. j .ne. 0 ) then
                  xp             = x0(ix,iy,k)
                  yp             = z0(ix,iy,k)
                  call area_subxz(i,j,k,xu,zu,xp,yp,as1,as2,as3,as4)
                  a_tot = as1 + as2 + as3 + as4
                  dene0(ix,iy,k,n) = 
     .            ( as1 * deneu(i,j,k,n)   + as2 * deneu(i+1,j,k,n) +
     .              as4 * deneu(i,j+1,k,n) + as3 * deneu(i+1,j+1,k,n) ) 
     .            / a_tot
                endif
             enddo
           enddo
         enddo
         do ip = 1,ipsmooth 
           call smoothx(dene0(:,:,:,n)) 
           call smoothy(dene0(:,:,:,n)) 
           call smoothz(dene0(:,:,:,n)) 
         enddo 
         write(34) dene0(:,:,:,n)
       enddo

!      initialize glon0 

       do k = 1,nl 
         do j = 1,ny
           do i = 1,nx
             glon0(i,j,k) = 0.
           enddo
         enddo
       enddo

!      calculate interpolated value of glon

       do k = 1,nl
         do iy = 1,ny
           do ix = 1,nx
              i              = i0(ix,iy,k)
              j              = j0(ix,iy,k)
              if ( i .ne. 0 .and. j .ne. 0 ) then
                xp             = x0(ix,iy,k)
                yp             = z0(ix,iy,k)
                call area_subxz(i,j,k,xu,zu,xp,yp,as1,as2,as3,as4)
                a_tot = as1 + as2 + as3 + as4
                glon1 = glonu(i,j,k)
                glon2 = glonu(i+1,j,k)
                glon3 = glonu(i,j+1,k)
                glon4 = glonu(i+1,j+1,k)
                glont = glon1 + glon2 + glon3 + glon4
                if ( (glon1 .gt. 355. .or.
     .                glon2 .gt. 355. .or.
     .                glon3 .gt. 355. .or.
     .                glon4 .gt. 355.     ) .and.
     .                glont .lt. 1200.            ) then
                  glon2 = glon1
                  glon3 = glon1
                  glon4 = glon1
                endif
                glon0(ix,iy,k) = 
     .          (   as1 * glon1
     .            + as2 * glon2
     .            + as4 * glon3
     ,            + as3 * glon4  )
     .          / a_tot 
              endif
           enddo
         enddo
       enddo

!      define longitude 'outside' of data

       j00 = 1
       do k = 1,nl 
         do i = 1,nx
           do j = 1,ny
             if ( glon0(i,j,k) .ne. 0 ) j00 = j
             if ( glon0(i,j,k) .eq. 0 ) 
     .            glon0(i,j,k) = glon0(i,j00,k)
           enddo
         enddo
       enddo

!      fix end points

       do k = 1,nl
         do j = 1,ny
           glon0(1,j,k)  = glon0(2,j,k)
           glon0(nx,j,k) = glon0(nx-1,j,k)
         enddo
       enddo

!      write out interpolated longitude file
  
       write(32) glon0

       stop
       end

       
*******************************************
*******************************************

!            area_subxz

*******************************************
*******************************************


        subroutine area_subxz(i,j,k,x,z,x0,z0,as1,as2,as3,as4)

!       calculate areas of cell sides
!       break each quadrilateral side into
!       two triangles and use the formula: 
!           A = (1/2)|a x b|
!       where
!           a: vector from A to B
!           b: vector from A to C

        include "param_diag.inc"

        real x(nz,nf,nl),z(nz,nf,nl)

!       as1

        ax1 = 0.5 * ( x(i+1,j,k) + x(i+1,j+1,k) ) - x0
        az1 = 0.5 * ( z(i+1,j,k) + z(i+1,j+1,k) ) - z0

        bx1 = x(i+1,j+1,k) - x0
        bz1 = z(i+1,j+1,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i,j+1,k) + x(i+1,j+1,k) ) - x0
        az2 = 0.5 * ( z(i,j+1,k) + z(i+1,j+1,k) ) - z0

        bx2 = x(i+1,j+1,k) - x0
        bz2 = z(i+1,j+1,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as1 =  a1 + a2  

!       as2

        ax1 = 0.5 * ( x(i,j+1,k) + x(i+1,j+1,k) ) - x0
        az1 = 0.5 * ( z(i,j+1,k) + z(i+1,j+1,k) ) - z0

        bx1 = x(i,j+1,k) - x0
        bz1 = z(i,j+1,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i,j,k) + x(i,j+1,k) ) - x0
        az2 = 0.5 * ( z(i,j,k) + z(i,j+1,k) ) - z0

        bx2 = x(i,j+1,k) - x0
        bz2 = z(i,j+1,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as2 =  a1 + a2  

!       as3

        ax1 = 0.5 * ( x(i,j,k) + x(i,j+1,k) ) - x0
        az1 = 0.5 * ( z(i,j,k) + z(i,j+1,k) ) - z0

        bx1 = x(i,j,k) - x0
        bz1 = z(i,j,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i,j,k) + x(i+1,j,k) ) - x0
        az2 = 0.5 * ( z(i,j,k) + z(i+1,j,k) ) - z0

        bx2 = x(i,j,k) - x0
        bz2 = z(i,j,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as3 =  a1 + a2  

!       as4

        ax1 = 0.5 * ( x(i,j,k) + x(i+1,j,k) ) - x0
        az1 = 0.5 * ( z(i,j,k) + z(i+1,j,k) ) - z0

        bx1 = x(i+1,j,k) - x0
        bz1 = z(i+1,j,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i+1,j+1,k) + x(i+1,j,k) ) - x0
        az2 = 0.5 * ( z(i+1,j+1,k) + z(i+1,j,k) ) - z0

        bx2 = x(i+1,j,k) - x0
        bz2 = z(i+1,j,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as4 =  a1 + a2  

        return
        end

******************************************************** 
******************************************************** 
******************************************************** 
******************************************************** 
 
        subroutine smoothx(f)  
  
        include "param_diag.inc" 
 
        parameter ( nnxp2 = nx + 2, nnyp2 = ny + 2, nnzp2 = nl + 2 )  
        parameter ( nnxp1 = nx + 1, nnyp1 = ny + 1, nnzp1 = nl + 1 )  
        parameter ( nnx   = nx    , nny   = ny    , nnz   = nl     )  
          
        real f(nnx,nny,nnz),f0(nnx+2,nny+2,nnz+2)  
  
        u12 = 1.  
  
        do k = 1,nnz  
          do j = 1,nny  
            do i = 1,nnx  
              f0(i+1,j+1,k+1) = f(i,j,k)  
            enddo  
          enddo  
        enddo 
  
!       zero-gradient in x  
 
        do k = 2,nnzp1  
          do j = 2,nnyp1  
            f0(1,j,k)     = f0(2,j,k)  
            f0(nnx+2,j,k)  = f0(nnx+1,j,k)  
          enddo  
        enddo 
  
!       zero gradient in y  
  
        do k = 2,nnzp1         
          do i = 1,nnxp2  
            f0(i,1,k)     = f0(i,2,k)  
            f0(i,nnyp2,k) = f0(i,nnyp1,k)  
          enddo        
        enddo 
 
!       periodic in z 
  
        do j = 1,nnyp2 
          do i = 1,nnxp2 
            f0(i,j,1)      = f0(i,j,nnzp1)  
            f0(i,j,nnzp2)  = f0(i,j,3)  
          enddo 
        enddo  
  
!       sweep in x (1/2)  
        
        do k = 2,nnzp1 
          do j = 2,nnyp1  
            do i = 1,nnxp1  
              f0 (i,j,k) = 0.5 * ( f0(i,j,k) + u12*f0(i+1,j,k) )  
            enddo  
          enddo  
        enddo 
  
        do k = 2,nnzp1 
          do j = 2,nnyp1  
            do i = nnxp2,2,-1  
              f0 (i,j,k) = 0.5 * ( f0(i,j,k) + u12*f0(i-1,j,k) )  
            enddo  
          enddo  
        enddo 
  
!       now get f  
  
        do k = 1,nnz 
          do j = 1,nny  
            do i = 1,nnx  
              f(i,j,k)  = f0(i+1,j+1,k+1)  
            enddo  
          enddo  
        enddo 
  
        return   
        end  

******************************************************** 
******************************************************** 
******************************************************** 
 
 
        subroutine smoothy(f)  
 
        include "param_diag.inc" 
 
        parameter ( nnxp2 = nx + 2, nnyp2 = ny + 2, nnzp2 = nl + 2 )  
        parameter ( nnxp1 = nx + 1, nnyp1 = ny + 1, nnzp1 = nl + 1 )  
        parameter ( nnx   = nx    , nny   = ny    , nnz   = nl     )  
        
        real f(nnx,nny,nnz),f0(nnx+2,nny+2,nnz+2)  
   
        u12 = 1.  
 
        do k = 1,nnz 
          do j = 1,nny  
            do i = 1,nnx  
              f0(i+1,j+1,k+1) = f(i,j,k)  
            enddo  
          enddo  
        enddo 
       
!       zero-gradient in x  
 
        do k = 2,nnzp1  
          do j = 2,nnyp1  
            f0(1,j,k)     = f0(2,j,k)  
            f0(nnx+2,j,k)  = f0(nnx+1,j,k)  
          enddo  
        enddo 
 
!       zero gradient in y  
  
        do k = 2,nnzp1         
          do i = 1,nnxp2  
            f0(i,1,k)     = f0(i,2,k)  
            f0(i,nnyp2,k) = f0(i,nnyp1,k)  
          enddo        
        enddo 
!       periodic in z 
  
        do j = 1,nnyp2 
          do i = 1,nnxp2 
            f0(i,j,1)      = f0(i,j,nnzp1)  
            f0(i,j,nnzp2)  = f0(i,j,3)  
          enddo 
        enddo  
  
!       sweep in y (1/2)  
 
        do k = 2,nnzp2 
          do j = 1,nnyp1  
            do i = 2,nnxp2  
              f0 (i,j,k) = 0.5 * ( f0(i,j,k) + u12*f0(i,j+1,k) )  
            enddo 
          enddo  
        enddo  
 
        do k = 2,nnzp1 
          do j = nnyp2,2,-1  
            do i = 2,nnxp2  
              f0 (i,j,k) = 0.5 * ( f0(i,j,k) + u12*f0(i,j-1,k) )  
            enddo  
          enddo  
        enddo 
 
!      now get f  
 
        do k = 1,nnz 
          do j = 1,nny  
            do i = 1,nnx  
              f(i,j,k)  = f0(i+1,j+1,k+1)  
            enddo  
          enddo 
        enddo  
 
      return  
      end  
******************************************************** 
******************************************************** 
******************************************************** 
 
        subroutine smoothz(f)  
  
        include "param_diag.inc" 
 
        parameter ( nnxp2 = nx + 2, nnyp2 = ny + 2, nnzp2 = nl + 2 )  
        parameter ( nnxp1 = nx + 1, nnyp1 = ny + 1, nnzp1 = nl + 1 )  
        parameter ( nnx   = nx    , nny   = ny    , nnz   = nl     )  
          
        real f(nnx,nny,nnz),f0(nnx+2,nny+2,nnz+2)  
  
        u12 = 1.  
  
        do k = 1,nnz  
          do j = 1,nny  
            do i = 1,nnx  
              f0(i+1,j+1,k+1) = f(i,j,k)  
            enddo  
          enddo  
        enddo 
  
!       zero-gradient in x  
 
        do k = 2,nnzp1  
          do j = 2,nnyp1  
            f0(1,j,k)     = f0(2,j,k)  
            f0(nnx+2,j,k)  = f0(nnx+1,j,k)  
          enddo  
        enddo 
  
!       zero gradient in y  
  
        do k = 2,nnzp1         
          do i = 1,nnxp2  
            f0(i,1,k)     = f0(i,2,k)  
            f0(i,nnyp2,k) = f0(i,nnyp1,k)  
          enddo        
        enddo 
!       periodic in z 
  
        do j = 1,nnyp2 
          do i = 1,nnxp2 
            f0(i,j,1)      = f0(i,j,nnzp1)  
            f0(i,j,nnzp2)  = f0(i,j,3)  
          enddo 
        enddo  
  
!       sweep in z (1/2)  
        
        do k = 1,nnzp1 
          do j = 2,nnyp1  
            do i = 2,nnxp1  
              f0 (i,j,k) = 0.5 * ( f0(i,j,k) + u12*f0(i,j,k+1) )  
            enddo  
          enddo  
        enddo 
  
        do k = nnzp2,2,-1 
          do j = 2,nnyp1  
            do i =2, nnxp1  
              f0 (i,j,k) = 0.5 * ( f0(i,j,k) + u12*f0(i,j,k-1) )  
            enddo  
          enddo  
        enddo 
  
!       now get f  
  
        do k = 1,nnz 
          do j = 1,nny  
            do i = 1,nnx  
              f(i,j,k)  = f0(i+1,j+1,k+1)  
            enddo  
          enddo  
        enddo 
  
        return   
        end  
 
