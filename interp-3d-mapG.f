
!
!      INTERP-3D-MAP.F
!
!      determine interior points for sami3 interpolations

       include "param_diag.inc"

!      arrays of 'real' data

       real glatu(nz,nf,nl),zaltu(nz,nf,nl)
       real xu(nz,nf,nl),zu(nz,nf,nl)
       real zualtmax(nf,nl)    !gj

!      arrays of interpolated data

       real    x0(nx,ny,nl),z0(nx,ny,nl)
       integer i0(nx,ny,nl),j0(nx,ny,nl)

       real glat0(nx,ny,nl),glon0(nx,ny,nl),zalt0(nx,ny,nl)

!      open grid data files

       open(unit=10,file='glatu.dat', form='unformatted')
       open(unit=11,file='zaltu.dat', form='unformatted')

!      open map files

       open(unit=34,file='map.dat', form='unformatted')

!      open interpolated grid files

       open(unit=35,file='glat0.dat', form='unformatted')
       open(unit=36,file='zalt0.dat', form='unformatted')
       open(unit=37,file='glon0.dat', form='unformatted')

       open(unit=45,file='x0.dat', form='unformatted')
       open(unit=46,file='z0.dat', form='unformatted')
       open(unit=55,file='xu.dat', form='unformatted')
       open(unit=56,file='zu.dat', form='unformatted')

!      some parameters

       re   = 6370.
       pi   = 4. * atan(1.)
       dtor = pi / 180.
       gamy = 3.

!      read in 'real' data

       read(10) glatu
       read(11) zaltu

!      get min and max of glatu and zaltu

       print *,'getting min/max glatu/zaltu - setting grid0'

       do k = 1,nl

         glatmin =   90.
         glatmax =  -90.
         zaltmin =   85.
         zaltmax =    0.

!      zaltmin is maximum altitude of lowest field line 
 
         j = 1 
         do i = 1,nz 
           if ( zaltu(i,j,k) .ge. zaltmin ) zaltmin = zaltu(i,j,k) 
         enddo 

         do j = 1,nf
           do i = 1,nz
             if ( glatu(i,j,k) .le. glatmin ) glatmin = glatu(i,j,k)
             if ( glatu(i,j,k) .ge. glatmax ) glatmax = glatu(i,j,k)
             if ( zaltu(i,j,k) .ge. zaltmax ) zaltmax = zaltu(i,j,k)
           enddo
           zualtmax(j,k) = zaltmax
         enddo

!        set up grid for glat0 (uniform)

         delg0  = ( glatmax - glatmin ) / float(nx-1)

         do iy = 1,ny 
           do ix = 1,nx
             glat0(ix,iy,k) = glatmin + ( ix - 1 ) * delg0
           enddo
         enddo

!        set up grid for glon0 (uniform)

         delg0  = 360. / float(nl-1)

         do iy = 1,ny 
           do ix = 1,nx
             glon0(ix,iy,k) = ( k - 1 ) * delg0
           enddo
         enddo

!        set up grid for zalt0 (non-uniform)

         zaltmax0 = 2000.

         do iy = 1,ny
           do ix = 1,nx
             iy0  = ny + 1 - iy
             dy   = float(iy-1) / float(ny-1)
             y2   = dy * sinh ( gamy )
             y3   = alog ( y2 + sqrt ( y2**2 + 1. ) ) / gamy
             zalt0(ix,iy0,k) = zaltmin + ( zaltmax0 - zaltmin ) 
     .                                           * ( 1. - y3 )
           enddo
         enddo

       enddo

!      obtain xu and zu 
!      (cartesian coordinates of glatu and zaltu)

       print *,'setting up xu/zu'

       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
             xu(i,j,k) = ( zaltu(i,j,k) + re ) * 
     .                 cos ( glatu(i,j,k) * dtor )
             zu(i,j,k) = ( zaltu(i,j,k) + re ) * 
     .                 sin ( glatu(i,j,k) * dtor )
           enddo
         enddo
       enddo

       write(55) xu
       write(56) zu

!        obtain x0 and z0 
!        (cartesian coordinates of glat0 and zalt0, i.e.,
!        the interpolated grid)

       print *,'setting up glat0/zalt0'

       do k = 1,nl
         do j = 1,ny
           do i = 1,nx
             x0(i,j,k) = ( zalt0(i,j,k) + re ) 
     .                   * cos ( glat0(i,j,k) * dtor )
             z0(i,j,k) = ( zalt0(i,j,k) + re ) 
     .                    * sin ( glat0(i,j,k) * dtor )
           enddo
         enddo
       enddo

       write(45) x0
       write(46) z0

!      determine interior point 
!      using cross-products from opposing vertices
!      (determines where interpolated grid data is
!      in terms of the 'read' grid data)

       print *,'getting interior points'

       do k = 1,nl
         print *,'k = ',k
         do iy = 1,ny
! limit the number of field lines searched.
! if zu is greater than any point on the "real" field line, we don't 
!  need to search it  -- gj
             jtmp = 1
!  Remember zalt0 is constant in ix
             do while(zalt0(1,iy,k) .gt. zualtmax(jtmp,k) .and. 
     .            jtmp .le. nf-1)
                jtmp = jtmp +1
             enddo
             j00 = max(jtmp-1,1)
             istart = 1
           do ix = 1,nx
             i0(ix,iy,k) = 0
             j0(ix,iy,k) = 0
             ifind   = 0

             do j = j00,nf-1   !-- gj
               if ( ifind .eq. 0 ) then
                 do i = istart,nz-1
                   if ( ifind .eq. 0 ) then
                     ax  = xu(i+1,j,k) - xu(i,j,k)
                     ay  = zu(i+1,j,k) - zu(i,j,k)
                     bx  = xu(i,j+1,k) - xu(i,j,k)
                     by  = zu(i,j+1,k) - zu(i,j,k)
                     cx  = x0(ix,iy,k)   - xu(i,j,k)
                     cy  = z0(ix,iy,k)   - zu(i,j,k)
                     axc = ax * cy - ay * cx
                     bxc = bx * cy - by * cx
                     if ( axc * bxc .le. 0 .and. 
     .                    axc .le. 0       .and.
     .                    ifind .eq. 0           ) then
                       ax  = xu(i+1,j,k) - xu(i+1,j+1,k)
                       ay  = zu(i+1,j,k) - zu(i+1,j+1,k)
                       bx  = xu(i,j+1,k) - xu(i+1,j+1,k)
                       by  = zu(i,j+1,k) - zu(i+1,j+1,k)
                       cx  = x0(ix,iy,k)   - xu(i+1,j+1,k)
                       cy  = z0(ix,iy,k)   - zu(i+1,j+1,k)
                       axc = ax * cy - ay * cx
                       bxc = bx * cy - by * cx
                       if ( axc * bxc .le. 0 .and.
     .                      axc .ge. 0             ) then
                         i0(ix,iy,k) = i
                         j0(ix,iy,k) = j
                         ifind   = 1
                         istart = max(1,i-2)
                       endif
                     endif
                   endif
                 enddo
               endif
             enddo
           enddo
         enddo
       enddo     

       write(34) i0,j0
       write(35) glat0
       write(36) zalt0
!       write(37) glon0

       end

       
