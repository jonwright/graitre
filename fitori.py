from WWMpeakassign import WWMcrystal, getrxryrz, WWMdiffractometer, \
    rotmatx, rotmaty, rotmatz, anglevecs2D
from ImageD11 import grain

import numpy as np
import pylab
import sys


energy = float(sys.argv[1])
pks  = sys.argv[2] 
ubi = sys.argv[3]
tilt = float(sys.argv[4])
TOLERANCE=float(sys.argv[5])

c=WWMcrystal()
c.set_wvln( 12.3985/ energy )
c.dmin = 0.6
u0 = grain.read_grain_file(ubi)[0].u # np.eye(3)
#print u0
rx, ry, rz = getrxryrz( u0 )
c.set_orientation(u0)
c.generate_hkls()
c.generate_gve()
d = WWMdiffractometer(c, pks)
d.axistilt = tilt
d.computek()
d.computeomegas()
print d.omega.shape


def getom(d, hkl, fri):
    c = d.crystal
    ub = np.dot( c.umatrix, c.bmatrix )
    g = np.dot( ub, hkl )
    #print "gs",g.shape,"ub",c.ub,"hkl",hkl
    ### begin computek ###
    kz = g[2]
    modg2 = (g*g).sum(axis=0)
    tilt = d.axistilt
    num = - modg2 - 2*tilt*kz
    den = 2*np.sqrt( 1 - tilt*tilt )
    kx = num/den
    arg = modg2 - kx*kx - kz*kz
    mask = arg < 0
    ky  = np.sqrt( arg ) * fri
    ### end computek ###
    omega = anglevecs2D( g[0], g[1], kx, ky )
    return np.degrees(omega)


h = d.pars.heights
order = np.argsort( h )
# there are 16 (111) peaks
print len(order)
print c.hkls.shape
#1/0
chosen = order[:]
oc = np.take(d.pars.centroids, chosen )
hc = np.take(d.pars.heights, chosen )
h2 = 1/(d.h*d.h).sum(axis=1)/3

import scipy.optimize
print "Fitting"
data = []

for i,o in enumerate(oc):
    # find closest peak
    dom = np.degrees( d.omega)
    dom = np.where( np.isnan(dom), -1e9, dom)
    j = np.argmin(abs(o-dom))
    if abs( o-dom[j] ) < TOLERANCE:
        data.append( ( d.h[j], d.fri[j], o ) )


SAVEFILE = False
def gof( p, args):
    global SAVEFILE
    rx, ry, rz, tilt, wvln = p

    d.axistilt = tilt
    d.crystal.set_wvln( wvln )
    umat = np.dot( np.dot( rotmatx(rx),rotmaty(ry)),rotmatz(rz))
    d.crystal.set_orientation( umat )
    s = []
    if SAVEFILE:
        sf = open("fitori.dat","w")
    for h,fri,obs in args:
        ocalc = getom( d, h, fri )
        s.append(obs-ocalc)
        if SAVEFILE:
            sf.write( "%-4d %-4d %-4d %-2d "%(h[0],h[1],h[2],fri))
            sf.write( "%.6f %.6f\n"%(obs,ocalc))
    return np.array(s)
def vc(c):
    for i in range(c.shape[0]):
        fac = 1/np.sqrt(c[i,i])
        c[i,:]*=fac
        c[:,i]*=fac
    return c

x0 = rx,ry,rz, d.axistilt,  12.3985/ energy
def f2(p, *args):
    g = gof( p, args )
    return np.dot( g,g )
print data[0]
print "Before",x0
#ret = scipy.optimize.fmin( f2, x0, data )
ret = scipy.optimize.leastsq( gof, x0, data, full_output=1 )
x, cov_x, info, mesg, ier = ret
SAVEFILE=True
s_sq = (gof(x, data)**2).sum()/(len(data)-len(x0))
cov_x *= s_sq
if ret in [1,2,3,4]:
    print "OK",mesg
print "After :",
for i in range(len(x)):
    print x[i], np.sqrt(cov_x[i,i])
l = x[-1]
E = 12.3985/l
dl =  np.sqrt(cov_x[-1,-1])
dE = E*dl/l
print "%.6f %.6f"%(E,dE)
print "Variance, covariance"
print vc(cov_x)
open("fitted.ubi", "w").write(("%.12f  %.12f  %.12f\n"*3)%tuple(d.crystal.umatrix.T.ravel()))
print "about to plot"
pylab.plot( [p[-1] for p in data], gof( x, data), "+" )
pylab.show()

print "end"
