module misc_mod

  use parameter_mod

!     diagnostic variables

    real :: u1(nz,nf,nl),u2(nz,nf,nl),u3(nz,nf,nl), &
            u4(nz,nf,nl),u5(nz,nf,nl)

    real :: ppt(nzp1,nfp1,nlt)
    real*8 :: blonp0t(nlt+3)

    real :: deni_mnp(nz,nion),ti_mnp(nz,nion),te_mnp(nz)

end module misc_mod
