
module photo_production_mod

  use parameter_mod

!     zenith data 
 
    real :: cosbdec(nz,nf,nl),sinbdec(nz,nf,nl),cx(nz,nf,nl), & 
            coschicrit(nz,nf,nl) 
 
!     photodeposition rates 
!     used 3 (absorption) and 7 (nion) explicitly 
!     used 4 (number of angles in nighttime deposition) 
 
    real :: sigabsdt(linesuv,3),flux(linesuv),sigidt(linesuv,7) 
    real :: sigint(linesnt,7),fluxnt(nz,nf,nl,91,linesnt) 
    real :: thetant(linesnt,4),zaltnt(linesnt,2) 
 
end module photo_production_mod
