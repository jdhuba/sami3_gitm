!     namelist data

module namelist_mod

  use parameter_mod

    logical :: hall,restart
    logical :: lmadala,lcr,lvs,lweimer,lhwm93,lhwm14

    integer :: psmooth,nion1,nion2
    integer :: maxstep,mmass

    real :: snn(7)
    real :: hrmax, dthr, hrpr, dt0, &
            rmin, altmin, fbar, f10p7, ap, &                            
            year, day, hrinit, tvn0, tvexb0, ver, veh, vw,&
            gamss, alt_crit, cqe, plat, plon,alt_crit_avg
    real :: storm_ti, storm_tf, vexb_max, &
            decay_time, pcrit, anu_drag0, &
            blat_max4, stn,denmin,blat_min

end module namelist_mod
