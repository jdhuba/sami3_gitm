
module gitm_mod

  integer,parameter ::  nLons = 90
  integer,parameter ::  nLats = 90
  integer,parameter ::  nAlts = 50
  real,parameter    ::  dt_gitm = 0.25


        real :: altitudei(nLons,nLats,nAlts)
        real :: uni(nLons,nLats,nAlts)
        real :: vni(nLons,nLats,nAlts)
        real :: temperaturei(nLons,nLats,nAlts)
        real :: n2i(nLons,nLats,nAlts)
        real :: o2i(nLons,nLats,nAlts)
        real :: oi(nLons,nLats,nAlts)
        real :: noi(nLons,nLats,nAlts)
        real :: n4si(nLons,nLats,nAlts)
        real :: hi(nLons,nLats,nAlts)
        real :: lat(nLats)
        real :: lon(nLons)

        real :: altitudef(nLons,nLats,nAlts)
        real :: unf(nLons,nLats,nAlts)
        real :: vnf(nLons,nLats,nAlts)
        real :: temperaturef(nLons,nLats,nAlts)
        real :: n2f(nLons,nLats,nAlts)
        real :: o2f(nLons,nLats,nAlts)
        real :: of(nLons,nLats,nAlts)
        real :: nof(nLons,nLats,nAlts)
        real :: n4sf(nLons,nLats,nAlts)
        real :: hf(nLons,nLats,nAlts)


end module gitm_mod
