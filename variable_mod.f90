!******************************************

module variable_mod

    use parameter_mod

    real :: deni(nz,nf,nl,nion),ne(nz,nf,nl),denn(nz,nf,nl,nneut)
    real :: denni(nz,nf,nl,nneut),dennf(nz,nf,nl,nneut)
    real :: vsi(nz,nf,nl,nion),vsid(nz,nf,nl,nion)
    real :: sumvsi(nz,nf,nl,nion),vsic(nz,nf,nl,nion)
    real :: te(nz,nf,nl),ti(nz,nf,nl,nion),tn(nz,nf,nl)
    real :: tni(nz,nf,nl),tnf(nz,nf,nl)
    real :: u(nz,nf,nl),v(nz,nf,nl),vpi(nz,nf,nl),w(nz,nf,nl)
    real :: ui(nz,nf,nl),vi(nz,nf,nl),wi(nz,nf,nl)
    real :: uf(nz,nf,nl),vf(nz,nf,nl),wf(nz,nf,nl)
    real :: vnq(nz,nf,nl),vnp(nz,nf,nl),vnphi(nz,nf,nl)

end module variable_mod
