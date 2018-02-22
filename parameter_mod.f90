
module parameter_mod

!      total number of workers
    integer,parameter ::  numwork = 6

!      number of altitudes (e.g., field lines)

    integer,parameter :: nf   = 124
    integer,parameter :: nfp1 = nf + 1
    integer,parameter :: nfm1 = nf - 1
    integer,parameter :: nfm2 = nf - 2 

!      number of grid cells along field line
!      nseg needs to be even
         
    integer,parameter :: nz0 = 200, nextra = 2, nseg = 2, &
                         nze = nseg * nextra, nz = nz0+nze, &
                         nzp1 = nz + 1, &
                         nzm1 = nz - 1  

!      number of grid cells in longitudinal direction per worker

    integer,parameter ::  nl   = 18, &
                          nlp1 = nl + 1, &
                          nlm1 = nl - 1  

!      Total number of grid cells in longitudinal direction

    integer,parameter ::  nlt   = numwork*(nl - 2), &
                          nltp1 = nlt + 1, &
                          nltm1 = nlt - 1  

!      grid for madala potential solver

    integer,parameter ::  nnx =  nlt + 1, nny = nf - 1 
    integer,parameter ::  nyextra = 5, nnyt = nny + nyextra 

!      ion densities

    integer,parameter :: nion  = 7    ! number of ions
    integer,parameter :: pthp  = 1    ! h+
    integer,parameter :: pthep = 5    ! he+
    integer,parameter :: ptnp  = 7    ! n+
    integer,parameter :: ptop  = 2    ! o+
    integer,parameter :: ptn2p = 6    ! n2+
    integer,parameter :: ptnop = 3    ! no+
    integer,parameter :: pto2p = 4    ! o2+

!      neutrals

    integer,parameter :: nneut = 7    ! number of neutrals
    integer,parameter :: pth   = 1    ! h
    integer,parameter :: pthe  = 5    ! he
    integer,parameter :: ptn   = 7    ! n
    integer,parameter :: pto   = 2    ! o
    integer,parameter :: ptn2  = 6    ! n2
    integer,parameter :: ptno  = 3    ! no
    integer,parameter :: pto2  = 4    ! o2

!      number of chemical reactions

    integer,parameter :: nchem = 21 

!      various constants

!      ftnchek is giving some meaningless errors about precision,
!      but i am going to lower the precision of some of these
!      variables to keep down the error messages


    real,parameter :: pie    = 3.1415927   
    real,parameter :: po180  = 1.745329e-02 
    real,parameter :: rtod   = 57.295780    
    real,parameter :: tm18   = 1.e-18       
    real,parameter :: charge = 4.8032e-10   

    real,parameter :: spday  = 86400., sphr = 3600.   
    real,parameter :: sol    = 3.e10        

    real,parameter :: gzero  = 980.665, re = 6370., bmag = 0.31 

    real,parameter :: bolt   = 1.38044e-16 
    real,parameter :: amu    = 1.67252e-24 
    real,parameter :: evtok  = 1.1604e+04  

!       real,parameter :: linesuv = 105, linesnt = 4   ! fism

    integer,parameter :: linesuv = 37, linesnt = 4     ! euvac

    real,parameter :: dayve = 80., sidyr = 365.4, solinc = 23.5 

!      these are for the error function

    real,parameter :: pas =   .3275911, z1 =  .2548295, &
                      z2  = - .284496 , z3 = 1.421413, &
                      z4  = -1.453152 , z5 = 1.0614054  

end module parameter_mod
