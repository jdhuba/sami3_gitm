module grid_mod

  use parameter_mod

    real :: alts(nz,nf,nl),grs(nz,nf,nl),glats(nz,nf,nl),glons(nz,nf,nl)
    real :: glons0(nz,nf,nl)
    real :: bms(nz,nf,nl),gs(nz,nl),ps(nz,nf,nl),gp(nz,nf,nl), &
            blats(nz,nf,nl),blons(nz,nf,nl),cfs(nz,nl)
    real :: gsthetax(nz,nf,nl),gsthetay(nz,nf,nl),gsthetaz(nz,nf,nl)
    real :: gsphix(nz,nf,nl),gsphiy(nz,nf,nl),gsphiz(nz,nf,nl)
    real :: gsrx(nz,nf,nl),gsry(nz,nf,nl),gsrz(nz,nf,nl)
    real :: ds(nz,nf,nl),d2s(nz,nf,nl),d22s(nz,nf,nl)
    real :: dels(nz,nf,nl)

    real :: xnorms(nzp1,nf,nl),ynorms(nzp1,nf,nl),znorms(nzp1,nf,nl)
    real :: xnormp(nz,nfp1,nl),ynormp(nz,nfp1,nl),znormp(nz,nfp1,nl)
    real :: xnormh(nz,nf,nlp1),ynormh(nz,nf,nlp1),znormh(nz,nf,nlp1)
    real :: xrg(nz,nf,nl),xthg(nz,nf,nl),xphig(nz,nf,nl)
    real :: qs(nz,nf,nl),brs(nz,nf,nl)

    real :: delhp(nz,nfp1,nl)
    real :: delph(nz,nf,nlp1)
    real :: delps(nzp1,nf,nl)
    real :: delhs(nzp1,nf,nl)

    real :: bdirhx(nz,nf,nlp1),bdirhy(nz,nf,nlp1),bdirhz(nz,nf,nlp1)
    real :: bdirsx(nz,nf,nl),bdirsy(nz,nf,nl),bdirsz(nz,nf,nl)
    real :: bdirpx(nz,nfp1,nl),bdirpy(nz,nfp1,nl),bdirpz(nz,nfp1,nl)
    real :: blatss(nfp1,nlp1),blonss(nfp1,nlp1),bradss(nfp1,nlp1), &
            blatsp(nf,nlp1),blonsp(nf,nlp1),bradsp(nf,nlp1), &
            blatsh(nfp1,nl),blonsh(nfp1,nl),bradsh(nfp1,nl)

    real :: vol(nz,nf,nl)
    real :: areap(nz,nfp1,nl),areas(nzp1,nf,nl),areah(nz,nf,nlp1)
    real :: vpnx(nzp1,nfp1,nlp1),vpny(nzp1,nfp1,nlp1),vpnz(nzp1,nfp1,nlp1)
    real :: vpsnx(nz,nf,nl),vpsny(nz,nf,nl),vpsnz(nz,nf,nl)
    real :: vhsnx(nz,nf,nl),vhsny(nz,nf,nl),vhsnz(nz,nf,nl)
    real :: vhnx(nzp1,nfp1,nlp1),vhny(nzp1,nfp1,nlp1),vhnz(nzp1,nfp1,nlp1)
    real :: xdels(nz,nfp1,nlp1),xdelp(nzp1,nf,nlp1),xdelh(nzp1,nfp1,nl)

    real :: blatp(nzp1,nfp1,nlp1),blatpt(nzp1,nfp1,nlt)
    real :: blonpt(nzp1,nfp1,nlt)
    real :: pp(nzp1,nfp1,nlp1),qp(nzp1,nfp1,nlp1),blonp(nzp1,nfp1,nlp1)
    real :: baltp(nzp1,nfp1,nlp1),brp(nzp1,nfp1,nlp1)
    real :: bmpf(nz,nfp1,nl),bmhf(nz,nf,nlp1),bmsf(nzp1,nf,nl)


end module grid_mod
