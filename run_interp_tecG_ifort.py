#!/usr/bin/python

import os,sys,string

try:
    rfile  = open('param_diag_blank.inc','r')
    wfile  = open('param_diag.inc','w')
    file   = rfile.read()
    print 'Number of time steps: '
    nt     = raw_input()
    line   = ('       parameter ( nt = %s ) \n' % nt)
    filen  = line+file
    wfile.write(filen)
    wfile.close()

    os.system( 'ifort -O3 -mcmodel=large -shared-intel -save interp-3d-mapG.f  -o interp-3d-mapG.x' )
    os.system( 'ifort -O3 -mcmodel=large -shared-intel -save interp-3d-dene_nG.f  -o interp-3d-dene_nG.x' )
    os.system( 'ifort -O3 -mcmodel=large -shared-intel -save tec_nmf2_hmf2-3dG.f  -o tec_nmf2_hmf2-3dG.x' )

    os.system( 'interp-3d-mapG.x' )    
    os.system( 'interp-3d-dene_nG.x' )    
    os.system( 'tec_nmf2_hmf2-3dG.x' )

except IOError:
    print 'File param_diag_blank.inc does not exist'
