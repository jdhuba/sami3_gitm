
module conductance_mod

  use parameter_mod

    real :: hipcpt(nf,nlt),hipcphit(nf,nlt),hihcmt(nf,nlt) 
    real :: hipct(nf,nlt),hihct(nf,nlt),hidvt(nf,nlt) 
    real :: hidpvt(nf,nlt),hidphivt(nf,nlt) 
    real :: hidpgt(nf,nlt),hidphigt(nf,nlt) 

    real :: hipc(nf,nl),hihc(nf,nl) 
    real :: hipcp(nf,nl),hipcphi(nf,nl),hihcm(nf,nl) 
    real :: hidv(nf,nl) 
    real :: hidpv(nf,nl),hidphiv(nf,nl) 
    real :: hidpg(nf,nl),hidphig(nf,nl) 
    real :: sigmap(nz,nf,nl),sigmah(nz,nf,nl) 
    real :: sigmapic(nz,nf,nl),sigmahic(nz,nf,nl) 


end module conductance_mod
