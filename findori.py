from WWMpeakassign import WWMcrystal, getrxryrz, WWMdiffractometer, \
    rotmatx, rotmaty, rotmatz
from ImageD11 import grain
import numpy as np
import sys

energy = 78.395
ubi  = "fitted.ubi"
pks  = sys.argv[1] # "wwm_setup_20160916_180a_.spc_pks.json"
NPKS = int(sys.argv[2])
c=WWMcrystal()
c.set_wvln( 12.3985/ energy )
c.dmin =5.43/np.sqrt(8.1)
#c.dmin = 
#gr = grain.read_grain_file( ubi )[0]
u0 = np.eye(3)
rx, ry, rz = getrxryrz( u0 )

c.set_orientation(u0)
c.generate_hkls()

c.generate_gve(  )

d = WWMdiffractometer(c, pks)
d.computek()
d.computeomegas()

h = d.pars.heights
order = np.argsort( h )
# there are 16 (111) peaks
chosen = order[-NPKS:]
oc = np.take(d.pars.centroids, chosen ) #- d.pars.zeroAngle
hc = np.take(d.pars.heights, chosen )

h2 = 1/(d.h*d.h).sum(axis=1)/3
#x,y =  np.load("www_plot.npy").T
x,y =  np.loadtxt("www_plot.dat").T

        

import pyqtgraph
pw = pyqtgraph.plot( title="My plot")
pw.plot( 180*d.omega/np.pi, h2, pen=None, symbol="+")
pw.plot( x, y )
pw.plot( oc, hc, pen=None, symbol="o")

import scipy.optimize
print "Fitting"


TRACK = 0
def func( x , *args):
    global TRACK
    rx, ry, rz = x
    #print "\tCalled with ",rx,ry,rz
    oc, d = args
    umat = np.dot( np.dot(  rotmatx(rx),rotmaty(ry)),rotmatz(rz))
    d.crystal.set_orientation( umat )
#    d.crystal.set_wvln( w )
    d.crystal.generate_gve()
    d.computek()
    d.computeomegas()
    dom = np.degrees( d.omega)
    esum = 0
#    print 'oc',oc
#    print 'dom',dom
    f = []
 #   print dom
    for i in range(len(oc)):
        err = (oc[i] - dom)
        j = (err*err).argmin()
        f.append( err[j] )
#        print i,oc[i],dom[j]
#    1/0    
    TRACK+=1
    return f

scormin=1e9
pmin=None
for i in range(1000): # if 1:
    if i>10 and np.random.random(1)<0.5:
        x0 = pmin + (np.random.random(3)-0.5)*np.pi/90
        before = scormin
    else:
        x0 = [0.,0.,0.]+(np.random.random(3)-0.5)*np.pi/2
        before = 0
    args = oc, d
    print "Start",x0,
    ret = scipy.optimize.leastsq( func, x0, args )[0]
    #def f2(x, *args):
    #    f = func(x, *args)
    #    return np.dot( f, f ).sum()
    #ret = scipy.optimize.fmin( f2, x0, args )
    print "ret",ret,
    fscor = func( ret, *args )
    scor = np.dot(fscor,fscor).sum()
    print "scor",scor,before
    if scor<scormin:
        scormin = scor
        pmin = ret
        





print
fscor = func( pmin, *args )    
scor = np.dot(fscor,fscor).sum()
print "BEST",pmin,"scores",scor
print fscor
print 12.3985/energy, energy
pw.plot( 180*d.omega/np.pi, h2, pen=None, symbol="x")
pw.plot( 180*d.omega/np.pi+180, h2+0.01, pen=None, symbol="+")

#for i,j in pairs:
#    pw.plot( [ 180*d.omega[j]/np.pi, oc[i]] ,
#             [ h2[j], hc[i]] )

raw_input()

open("fitted.ubi", "w").write(("%.12f  %.12f  %.12f\n"*3)%tuple(d.crystal.umatrix.T.ravel()))
