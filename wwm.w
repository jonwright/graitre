% nuweb formatted latex document for ID11 processing.
% 
% Please stay within 80 characters
% 
% Copyright: Jon Wright, 2013, ESRF


\documentclass[11pt,notitlepage]{article}
\usepackage{amsmath}
\usepackage{graphicx}

\usepackage[margin=2cm,a4paper]{geometry}

\begin{document}

\title{WWM : Wafer Wavelength Monitor, collected notes}
\author{Jon Wright}
\date{Summer 2013}
\maketitle

\begin{abstract}
This is the complete documentation package for the wavelength wafer 
monitor project.
A silicon wafer is scanned in angle in the direct beam and the transmission 
is recorded.
Code for the data collection is included here.
Whenever diffraction occurs an extinction peak is observed, and the
angles depend on the crystal orientation and wavelength.
Data reduction and wavelength fitting code for such data are developed here.
\end{abstract}

\tableofcontents

\section{Introduction}

A fast scan has been developed to measure the extinction diffraction
of a silicon wafer in order to make a wavelength monitor for ID11.
Data are recorded using the ESRF MUSST card for the incident and 
transmitted beam intensity and these are written to an ascii spec file.

In order to extract the wavelength, crystal orientation and axis direction
from the data some processing is required. The following steps are identified:
\begin{itemize}
\item Normalisation of the scan data
\item Extraction of peak positions, widths and heights as a function of angle
\item Assignment of hkl indices to the peaks
\item Refinement of the geometry (3 orientation, one axis tilt, one wavelength)
\end{itemize}


\section{Some scattered notes on extinction diffraction}

We are concerned by the overall length of the diffraction vector and
the projection of this perpendicular to the rotation axis.
For a single peak we can observe both $g$ and $-g$ as they each cross the
two sides of the Ewald sphere (from the inside to the outside and the 
reverse).
From these four observations we will find peaks at $\phi$ and
$\phi+180$ when the rotation axis is perpendicular to the beam.
If there is a beam tilt then there is an offset which depends on the azimuthal
angle of the reflection.
The difference in angle between two Friedel pairs gives us a single 
observation to constrain the length and axis projection (two unknowns).
If the length is known from having the hkl indices, cell parameter
and wavelength then we can determine the projection on the axis.

If we tweak the wavelength we change the length of the scattering vector 
but keep the same direction cosines. 
This tells us the sign of the solution in the square root equation as the
peaks shift either to lower or higher angles. 

Example data from 65 and 65.1keV show this in action (figure please FIXME).
Peaks shift by their width (about 100 eV) either to positive or negative.
Looks like some Freidel pairs get closer together in angle, some further apart.
None of them look like staying where they are...

There is a flick up on a tail of the peak sometimes (figure please).
Currently this is interpreted as a dynamical diffraction effect.
For the beam going directly through the crystal there is the usual 
absorbtion along the path.
If the wafer is at a grazing angle the diffracted beam can exit
on the opposite surface with a much smaller path length.
When this diffracted beam diffracts again it will have found a "short cut"
through the crystal.
To compute the peakshape for this properly we should look into what is 
going on inside the Bormann fan. 
(An interim solution is to avoid using such peaks).

There is also the anomalous transmission effect, when the standing 
wave inside the crystal avoids the atomic positions (see  Batterman and
Cole, 1964). 
To quantify these various effects we should compute the actual values 
of $\mu t$ as a function of energy and also figure out what are 
the Pendelossung periods etc.



\section{Axis tilt alignment procedure}

If a tilt is installed under the axis measure a pair of sharp peaks at 
0 and 180 degrees.
Plot the separation of these peaks (-180) versus the tilt angle.
When the axis is perpendicular to the beam the peaks will be 180 degrees 
apart.
Example scan data were collected, for 0,1,2,3 degrees these are:
\begin{verbatim}
lid112:/data/id11/crystallography/jon/june13/wwm_35kev
total 358740
-rw-r--r--  1 opid11 id11  1119343 Jul 30 19:10 junk.spc
-rw-r--r--  1 opid11 id11 66398182 Jul 30 15:32 wwm_35kev.spc
-rw-r--r--  1 opid11 id11 10885993 Jul 30 16:09 wwm_42kev_miome.spc
-rw-r--r--  1 opid11 id11 70560495 Jul 30 17:10 wwm_42kev_miome_hry.spc
-rw-r--r--  1 opid11 id11 70579160 Jul 30 17:31 wwm_42kev_miome_hry1.spc
-rw-r--r--  1 opid11 id11 70574330 Jul 30 17:48 wwm_42kev_miome_hry2.spc
-rw-r--r--  1 opid11 id11 70540335 Jul 30 19:04 wwm_42kev_miome_hry3.spc
-rw-r--r--  1 opid11 id11  5200492 Jul 30 16:47 wwm_42kev_miome_speed.spc

lid112:/data/id11/crystallography/jon/june13 % 
grep "#[OP]3" wwm_energies/wwm_hry_align_proc.spc/wwm_hry_align_proc.spc.spc 
#O3     tabx      taby    eurosp     cdtez   berger1     HrotY      s6vg      s6vo
#P3 3.051758e-07 0 600 -46.22 0  0   0.70000016 1.2903226e-07
#P3 3.051758e-07 0 600 -46.22 0  0   0.70000016 1.2903226e-07
#P3 3.051758e-07 0 600 -46.22 0  0.1 0.70000016 1.2903226e-07
#P3 3.051758e-07 0 600 -46.22 0  0.1 0.70000016 1.2903226e-07
#P3 3.051758e-07 0 600 -46.22 0 -0.1 0.70000016 1.2903226e-07
#P3 3.051758e-07 0 600 -46.22 0 -0.1 0.70000016 1.2903226e-07
#P3 3.051758e-07 0 600 -46.22 0 -0.063194444 0.70000016 1.2903226e-07
#P3 3.051758e-07 0 600 -46.22 0 -0.063194444 0.70000016 1.2903226e-07


lid112:/data/id11/crystallography/jon/june13 % more ~opid11/hry.dat 
0.0 0.02234
0.1 0.058066
-0.1 -0.013028

\end{verbatim}


\section{Thoughts on indexing crystals}

For the wavelength monitor we can just assume the crystal
orientation is known.
In practice one can read off the positions of the (111) reflections
as the strongest ones in the pattern anyway for a silicon wafer.

Nevertheless it is interesting to think about how to deal with an
unknown unit cell and to determine the orientation.

If the rotation is on a tilt stage and we can tilt the axis
with respect to the beam and the reflections shift according to their
azimuthal position. 
This should offer a way to get the (sin of) the azimuth angle.

If we tweak the wavelength the reflections generally shift 
either to higher or lower angles. 
This tells us whether the reflection crosses into or out of the
Ewald sphere as the rotation axis moves.
It also should tell us which direction to look in to find the 
Friedel pair.

The width of the peaks, at ID11, tells us the glancing angle
between the reflection and the thickness of the Ewald sphere 
for that reflection.
This depends on the bandpass and beam divergence as well as the
glancing angle. 

By combining several of the pieces of information 
one can presumably find orientation matrices in simple cases. 
Due to the symmetry of the experiment it is not clear whether
a solution is unique...


\section{Theory}

There is a nice derivation of the angle of diffraction in the 
rotation method in Milch and Minor (1974).
This is also something like the Bond method (ref). 
We reproduce the Milch and Minor story
with the notation changed to be more similar to fable/ImageD11.
The condition for diffraction is:

\[ \mathbf{| s_0 + k | = |s_0| } \]

Here $\mathbf{s_0}$ is the incident wave-vector and $\mathbf{k}$ is the 
crystal lattice vector at the moment when the diffraction occurs.
This equation just says the incident momentum is the same as the scattered.

\[ \mathbf{| s_0 + k |^2 = |s_0|^2 } \]
\[ \mathbf{ |s_0|^2 + 2s_0.k + |k|^2 - |s_0|^2 = 0}  \]
\[ \mathbf{  |k|^2 + 2s_0.k} = 0 \]

Once the hkl is assigned to the peak then $\mathbf{|k|^2 = |g|^2}$ is a 
constant which does not depend on rotations ($\mathbf{g}$ is the crystal
lattice vector prior to rotation to the diffraction condition).

We define a co-ordinate system where the rotation axis defines z and the 
beam is incoming roughly along x in the x-z plane.
This follows Milch and Minor (1974) but exchanges labels to approximately 
follow the fable convention (REF) with a key difference; 
\emph{the z axis is parallel to the rotation axis and units make 
$\mathbf{|s_0|}=1$}.
During the experiment the component of the vector along the rotation axis,
$g_z$, does not change, so that the rotated vector component 
$k_z$ is also constant.

Ideally, if the beam is perpendicular to the axis, then it will be along x.
Otherwise we have:

\[ \mathbf{s_0 = s_{0x} + s_{0z}}\]

...where since $ \mathbf{s_0} = 1 $ we can figure out that 
 $ s_{0x} = \sqrt{ 1 - s^2_{0z}}$.  
If we put these equations together we get:

\[\mathbf{  |k|^2 + 2(s_{0x} + s_{0z}).(k_x + k_y + k_z) = 0 } \]

Due to the choice of co-ordinate system and use of orthogonal axes 
only the xx, yy and zz terms from the dot product can be non-zero:

\[ \mathbf{ |k|^2 + 2 s_{0x}.k_x + 2 s_{0z}.k_z = 0 } \]
\[ \mathbf{ 2 s_{0x}.k_x = - |k|^2 - 2 s_{0z}.k_z } \]
\[ \mathbf{ k_x = - \frac{|k|^2 - 2 s_{0z}.k_z }{ 2 s_{0x}}} \]
\[ \mathbf{ k_x = - \frac{|k|^2 - 2 s_{0z}.k_z }{ 2 \sqrt{1-s^2_{0z}}}} \]
...via the choice of z along rotation axis then...:
\[ \mathbf{ k_x = - \frac{|g|^2 - 2 s_{0z}.g_z }{ 2 \sqrt{1-s^2_{0z}}}} \]

If the orientation matrix and beam direction are known we are able to 
compute the g-vector and ${s0z}$ and we can compute $k_x$.

Should someone take the time to carefully align the axis to be exactly
perpendicular to the beam then the $s_{0z}$ is zero and this reduces to:

\[ \mathbf{ k_x = -\frac{|g|^2}{2}} \]

Since the component along the axis is constant, $\mathbf{k_z=g_z}$, 
then we can get the
$\mathbf{k_y}$ component from $k_y=\pm \sqrt{|k|^2-k_x^2-k_z^2}$.

There are generally two positions of $\mathbf{k}$ where the diffraction 
can happen, corresponding to $\pm\mathbf{k_y}$.
We find these by solving the goniometer equation:

\[ k_x = g_x \cos{\phi} - g_y \sin{\phi} \]

Re-write this in terms of a some different quantities to make it possible
to solve for the angle. 
We use the length of the vector in the $xy$ plane, $g_{xy}$ and the angle
$\theta$:
\[g_{xy} = \sqrt{g_x^2+g_y^2};   \theta = \arctan{ g_y/g_x } \]
\[g_{x} = g_{xy}\cos{\theta}; g_y = g_{xy}\sin{\theta} \]
\[ k_x = g_{xy}\cos{\theta}\cos{\phi} - g_{xy}\sin{\theta}\sin{\phi} \]
\[ k_x = g_{xy}\cos{( \theta+\phi ) } \]
\[ \cos{(\theta+\phi)} = -\mathbf{ \frac{|g|^2 - 2 s_{0z}g_z }
     { 2  \sqrt{g_x^2+g_y^2} \sqrt{1-s_{0z}^2 }}} \]

Finally:

\[ \phi = \arccos{ \left[ -\frac{|g|^2 - 2 s_{0z}.g_z }
 { 2 g_{xy} \sqrt{1-s^2_{0z}}}\right]} - \arctan{ g_y/g_x } \]

There are two solutions for the inverse cosine for arguments
inside the range -1 to 1, 
no solutions outside the range and a single solution for $\pm 1$.

\subsection{Computing all this with derivatives}

We assume an initial orientation matrix is known and that we want to refine
orientation updates, the wavelength and axis tilt.
We write our crystal scattering vector as:

\[ \mathbf{c = U.B.h } \]

Going to the units for $g$ for this document means also multiplying 
by the wavelength, $\lambda$:

\[ \mathbf{g} = \lambda \mathbf{ U.B.h } \]

We need the derivative of $\mathbf{g}$ with respect to small orientation
changes, which we can choose to make on our $x,y,z$ axes.
After the small change, for example $r_z$ about $z$ we have:

\[ \mathbf{g} = (g_x\cos{r_z} + g_y\sin{r_z})\mathbf{i} +
   	        (-g_x\sin{r_z} + g_y\cos{r_z})\mathbf{j} +
		g_z\mathbf{k} \]
\[ \frac{d\mathbf{g}}{dr_z} = 
   (-g_x\sin{r_z} + g_y\cos{r_z})\mathbf{i} +
   (-g_x\cos{r_z} - g_y\sin{r_z})\mathbf{j} \]

We evaluate this at $r_z=0$ to obtain:

\[ \frac{d\mathbf{g}}{dr_z} = g_y.\mathbf{i} - g_x.\mathbf{j} + 0.\mathbf{k} \]

...which is just the vector $[g_y,-g_x,0]$. 
The same can be done for each of the three axes. 
These are just the cross products between the g-vector the three
basis vectors.

For the derivative with respect to the wavelength we have:

\[\frac{d \mathbf{g}}{d\lambda} = \mathbf{ U.B.h } = 
  \frac{\mathbf{g}}{\lambda} \]

For the derivatives of $|g|^2$, just $2\times$ a dot of the derivative 
with the original vector:

\[ \frac{ \partial{\mathbf{|g|^2}} }{ \partial x } = 
   \frac{ \partial{\mathbf{u.v}} }{ \partial x } = 
   \frac{ \partial{\mathbf{u^T v}} }{ \partial x } = 
   \frac{ \partial{\mathbf{u}} }{ \partial x }\mathbf{v} + 
   \frac{ \partial{\mathbf{v}} }{ \partial x }\mathbf{u} = 
   2 \frac{ \partial{\mathbf{g}} }{ \partial x }.\mathbf{g} 
 \]

For the derivatives of the inverse cosine and tangent terms we have:

\[ \frac{ \partial atan2(y,x) } { \partial x } = 
   \frac{ \partial \arctan(y/x) } { \partial x } = 
   - \frac{ y } { x^2 + y^2 } \]
\[ \frac{ \partial atan2(y,x) } { \partial y } = 
   \frac{ \partial \arctan(y/x) } { \partial y } = 
   \frac{ x } { x^2 + y^2 } \]
\[ \frac{ \partial \arccos(x) } { \partial x } = 
   -\frac{1}{ \sqrt{1-x^2} } \]

This should be enough mathematics to put together a fitting program.
The discontinuity at certain angles might be an issue. 
Since we work with a wafer we can design the program to place the 
problem at angles where there is no diffraction,
Our internal zero angle, and cut point should be close to the 
plane of the wafer to avoid this problem.


\subsection{ Experimental calibration, determination of $\mu$ }

In the ideal case the diodes and electrometers should be calibrated
to know the dark current and intensity when there is no crystal in 
the beam.
If this is not done we can try to figure out some values from the 
data.
In the test experiments scans were done where the rotation
axis was moved out of the beam, so we see the direct beam 
intensity and can measure it.
We also measured with the shutter closed to get the 
diode offset values.
The absolute X-ray flux can be determined from the diode
values via the method in Owen et al (2008).

Using the Beer-Lambert law the intensity transmitted will be:

\[    I = I_0 \exp{(-\mu l)} \]

...where $l$ is the path length through the wafer. 
The effective thickness of the wafer when it rotates is given by:

\[ l = \frac{t}{|\sin{(\phi+\phi_0)}|} \]

...where we have defined the zero angle to be where the wafer
is parallel to the beam.
If we take logs of the Beer-Lambert law and insert the equation
for $l$ and rearrange we get:

\[  |\sin{\phi+\phi_0}| (\log{I_0} - \log{I}) = \mu t   \]

Potential problems arise here if we do not have the true 
value for $I$, the intensity, but instead we measure some
number which is related to $I$ via a scale factor (gain) and
offset (dark current):

\[ I = s(V-V_d) \]

...where $s$ is a scale factor and $V, V_d$ are the output voltages
measured with and without the X-ray beam.


\[   \frac{I}{I_0} = \exp{(-\mu l)} =
   \frac{s(V-V_d)}{s_0(V_0-V_{0d})} \]

\[   \log{I_0} - \log{I} = \frac{\mu t}{|\sin{(\phi+\phi_0)}|} =
 \log{(s_0(V_0-V_{0d}))} - \log{(s(V-V_d))}  \]

\[ \frac{\mu t}{|\sin{(\phi+\phi_0)}|} = \log{(V_0/V)} + \log{(s_0/s)} + 
   	     \log{(1-V_{0d}/V_0)} - \log{(1-V_d/V))}  \]

Taking $log(1+x)=x$ when $x$ is small and writing $S=\log{(s_0/s)}$ gives:

\[ \frac{\mu t}{|\sin{(\phi+\phi_0)}|} ~= \log{(V_0/V)} + S - V_{0d}/V_0 - V_d/V  \]

\[ \frac{\mu t}{|\sin{(\phi+\phi_0)}|} ~= 
   \log{(V_0/V)} + S - \frac{ V_{0d}+V_d }{V_0V}  \]

We should check numerically, but it looks like the term in the dark currents
is likely to be much smaller than all the others.
So on a plot of $1/|sin{(\phi+\phi_0)}|$ we should have gradient $\mu t$ and 
intercept close to $S = \log{(s_0/s)}$.
If the $phi_0$ value is not correct there is a strong difference between
the $\mu t$ values for different parts of the scan and this is a 
non-linear effect. 

We need a little fitting program which takes the data and fits $\phi_0$ and
the scale factors.
Stepwise you can get a first guess of $\mu t$ from the gradients and adjust the
small angle regions to match by varying $\phi_0$.

It looks like a plot of the variation from the linear fit gives a 
straight line when plotted as:
\begin{verbatim}
plot(sign(tan(radians(msa)))*x, (y-np.polyval(p,x))/x,",")
\end{verbatim}
This should allow the zero angle to be determined.

Probably this just needs some least squares parameter fitting?

If we plot $\log{(V_0/V_k)}$ versus $1/|\sin(\phi+phi_0)|$ then
we should have a straight line with gradient $\mu t$ and 
intercept which picks up a lot of the crap above.

This trial has an "ang" definition which puts the zero at 90 degrees, anyway:

\begin{verbatim}
plot( 1/(abscosd(ang[::10])), 
      log( data[120400::10,3])-log(data[120400::10,2]),",")
ylim(0,1.5)
xlim(0,40)
title("wwm_42kev/wwm_42kev_0.0016.spc")
xlabel("1/|sin(phi+phi0)|")
ylabel("log(Vmonitor)-log(Vsignal)")
\end{verbatim}

This is plotted on figure~\ref{fig:muplot} on page~\pageref{fig:muplot}. 
It is clearly linear for a wide range of angles, up to the point where
the beam starts to overspill at the sides of the wafer.

If the $\phi_0$ value is not well determined we will get quite different
gradients in the different quadrants of 

\begin{figure}[tb]
\label{fig:muplot}
\includegraphics[width=\textwidth]{determine_mu}
\caption{ By measuring the gradient from this plot we can determine
the instrumental constants relating to diode efficiency and $\mu t$.}
\end{figure}

We note experimentally that fitting is mostly sensitive to the 
$\phi_0$ at angles near where the wafer is parallel to the beam, eg, 
near 0 and 180 degrees when $\sin{\phi}$ is also small.
We will assume that the offset angle $\phi_0$ that we are trying to 
correct is small.
In an iterative procedure this will eventually become true as the 
algorithm converges.
Expanding $\sin(\phi+\phi_0) = \sin{\phi}\cos{\phi_0} + \cos{\phi}\sin{\phi_0}$




\section{Fitting software}

The fitting code will be benchmarked against ImageD11 in the case of 
having no axis tilt. 
Following the easy maintenance of id31sum in comparison to anything
else the author has ever written, we try for a similar programming
style here.
That is to say all the modules are in the section in this nuweb 
document and the entire programs will be derived from this 
single file, probably as large single files themselves. 

In terms of overall structure we have a couple of approaches. 
We can correct the raw data and extract peak positions, then 
do a fit again the extracted positions.
This approach should be the first step.
Having gotten that working we can also do a fit against 
the raw data, Rietveld style, using the parameters determined.

This breaks up into several distinct steps.
First we take the raw data and find out the mu and wafer zero 
offset by fitting the absorbtion. 
Then we extract the peak positions, widths and heights.
Then we fit the positions and assign hkl indices.
Then we can go back into the data with a more complete hkl list.

@subsection{General coding stuff}

First the module imports. 
As this list grows so the distribution of the software is gradually
crippled.

@d imports
@{
import math, numpy as np, json, pprint
from ImageD11 import unitcell, gv_general, transform
from ImageD11 import connectedpixels
@}


Some little utility scripts for angle mods etc.

@d utilities
@{
def angmod(x): 
    """ Angle mod 360, radians """
    return np.arctan2( np.sin(x), np.cos(x))

def angmoddeg(x):
    """ Angle mod 360, degrees """
    return np.degrees(angmod(np.radians(x)))

def abscosdeg(x):
    """ Absolute cosine of angle in degrees """
    return np.abs(np.cos(np.radians(x)))

def abssindeg(x):
    """ Absolution sin of angle in degrees """
    return np.abs(np.sin(np.radians(x)))

@}


We will want some sort of description of a peak. 
This could be a python object, or the list of values to expect
in a C struct, or something else. 
For the moment we used connectedpixels from ImageD11 and this 
was set up to return arrays of centroids, areas, widths and 
heights.






Collect together all the little functions etc into one place for now:

@d wwmfunctions
@{
@< imports @>
@< utilities @>
@< generatehkls @>
@}

In the future we might want to split this into a no-dependency
version.

\subsection{Reading the experimental data}

The format is (unfortunately) defined in the wwm\_save spec
macro, which is not really ideal. 
In the near future we should save the diode dark offsets and 
full scale values in headers or a calibration file.


@d readspc
@{

def readspc(fname):
    """ minimalist spec file reading """
    scans = []
    data  = []
    coltitles = []
    nvals = -1
    for line in open(fname).xreadlines():
#        print line
        if line[0:2] == "#S":
            if len(data)>0:
                scans.append(np.array(data))
                data = []
            continue
        if line[0:2] == "#L":
            titles = line[2:].strip().split("  ") 
            nvals = len(titles)
    	    coltitles.append(titles)
            continue
        if line[0] == "#":
            continue
        vals = line.split()
        if len(vals)>0 and len(vals) == nvals:
            data.append( [float(v) for v in vals] )
    if len(data)>0:
        scans.append(np.array(data))
    return scans, coltitles

@}


\subsection{ Getting the signal from the raw data }

We need to define some experimental parameters.
We use a simple python dictionary and wrap this into
an object giving load and save capability and experiment description.

@d WWMpar
@{
wwmpars = {
    'zeroMonit' : 0.0,
    'zeroTrans' : 0.0,
    'zeroAngle'  : 0.0,
    'scaleMonit' : 1.0,
    'scaleTrans' : 1.0,
    'scaleAngle' : 160000.0,
    'chanMonit'  : 'CH6',
    'chanTrans'  : 'CH5',
    'chanAngle'    : 'CH1',
    }


class WWMpar(object):
    def __init__(self, filename=None, **pars):
        if filename is not None:
            self.load(filename)
        self.__dict__.update(pars)
    def save(self, filename):
        open(filename,"w").write( 
            json.dumps( self.__dict__,
            sort_keys=True,
            indent = 4,
            separators = (",",": " ))
        )
    def load(self, filename):
        d = json.loads(open(filename,"r").read())
        self.__dict__.update( d )
    def __str__(self):
        return str(self.__dict__)
@}

This little script checks we can load and save parameters
into a json file (ascii representation of a dict).

@O partest.py
@{
@<imports@>
@<WWMpar@>

if __name__=="__main__":
    p = WWMpar( **wwmpars )
    print str(p)
    p.save("testpars.json")
    p.load("testpars.json")
    print str(p)
@}

Assuming we have read in one or several spec files we need to 
fit them to extract the $\mu t$ and wafer zero angle.



\subsection{ Programs }



We generate idealised data with a (semi-) random orientation matrix
and compute the g-vectors and omega rotation angles. 
This little program does a simulation with a slightly tilted crystal.


This script computes the omega angles using the xfab/ImageD11 code.

@d generatehkls
@{
def generatehkls( a, wvln ):
    """
    Generate all of the hkls which can be reached 
    a is the cubic unit cell parameter
    """
    uc = unitcell.unitcell( [ a,a,a,90,90,90],"F" )
    uc.gethkls_xfab( 2.0/wvln, "Fd-3m" ) # out to 311
    return np.array([p[1] for p in uc.peaks])
@}


@o simulate_wwm.py
@{
@< wwmfunctions @>

if __name__=="__main__":
   a = 5.43094
   wvln = 0.25
   hkls = generatehkls( a, 2*a/np.sqrt(11.1)) # out to 311 for demo

   r3 = 1/math.sqrt(3)
   u = gv_general.rotation_axis([r3,r3,r3],2.0).matrix
   u = np.dot(u, gv_general.rotation_axis([0,0,1],45.0).matrix)
   ub = u/a     # wavelength in ImageD11 is supplied later
   gve = np.dot(ub, hkls.T )
   sol1, sol2, valid = gv_general.g_to_k( gve, wvln, axis=[0,0,-1] )
   tth, eta, omega = transform.uncompute_g_vectors(
   	  gve, wvln)
   omega = angmoddeg(omega)
   sol1 = angmoddeg(sol1)
   sol2 = angmoddeg(sol2)
   print "# U matrix"
   print "U = ",repr(u)
   print "wvln = ", wvln
   print "a =",a
   print "#  h   k   l   gx       gy       gz        om1      om2     ok  prec"
   for i in range(len(gve[0])):
        print ("%4d"*3)%tuple(hkls[i]),(" %8.5f"*3)%(tuple(gve[:,i])),
	if valid[i]: 
            print (" % 6.2f"*2)%(sol1[i],sol2[i])," valid",
        else:
            print "   NaN    NaN   invalid",
        dsame = abs(omega[0][i] - sol1[i]) + abs( omega[1][i]-sol2[i] )
        doppo = abs(omega[0][i] - sol2[i]) + abs( omega[1][i]-sol1[i] )
        print "%.1g"%(min(dsame,doppo))
@}

Numerically we could evaluate derivatives using the code above
for some parameters: ...


\subsection{ The user interface }

This is where all the time will go.

What is needed?

Ability to select a spec file.

Ability to decide on a scan and regions to use for dark and gain calibration
when wafer was offset.

Ability to find gain and offset by some other method (eg: noise in the monitor).

Ability to select scan(s) ready for peaksearching.

Ability to plot data and edit parameters.

Save peaks into a file.

Fit peak positions.

Globally we have a set of experimental parameters in the WWMexp object.
Separate from this we have a notion of a particular scan or series 
of scans giving rise to a dataset.

@d WWMdataset
@{

@<readspc@>
class WWMdataset(object):
    def __init__(self, fname, pars):
        self.fname = fname
        self.scans, self.coltitles = readspc( fname ) 
        self.pars = pars
    def getcol(self, scan, name):
        """ Accepts a single scan or a list of scans """
        if np.iterable(scan): 
            # Bounds check first
            i = self.coltitles[scan[0]].index(name)
            for s in scan:
                assert s >= 0 and s < len(self.scans), \
                    (s,scan, len(self.scans))
            assert i == self.coltitles[s].index(name)
            return np.concatenate( [self.scans[s][:,i] for s in scan] )
        else:
            if name in self.coltitles[scan]:
                i = self.coltitles[scan].index(name)
            else:
                print name
                print self.coltitles[scan]
                raise
            return self.scans[scan][:,i]
    def getMonit(self, scan):
        return self.getcol( scan, self.pars.chanMonit )
    def getTrans(self, scan):
        return self.getcol( scan, self.pars.chanTrans )
    def getAngle(self, scan):
        return self.getcol( scan, self.pars.chanAngle )/self.pars.scaleAngle
    def getSignal(self, scan):           
        t = (self.getTrans(scan) - self.pars.zeroTrans)/self.pars.scaleTrans
        m = (self.getMonit(scan) - self.pars.zeroMonit)/self.pars.scaleMonit
        return t/m


@}

\subsection{Dark and Gain Calibration}

We combine a set of experimental parameters with a scan to 
get to the calibration. 
Minimal inputs are the specfilename and scan numbers, we can also 
add the scan motor steps and angle channel.

@d WWMcalib
@{
@<UI@>

def WWMcalib( dataset, darkscan, floodscan ):
    dataset.pars.zeroMonit = dataset.getMonit( darkscan ).mean()
    dataset.pars.zeroTrans = dataset.getTrans( darkscan ).mean()
    fm = dataset.getMonit( floodscan ) - dataset.pars.zeroMonit
    ft = dataset.getTrans( floodscan ) - dataset.pars.zeroTrans
    print "Going to call UIgetrange"
    low, high = UIgetrange( "Range to use for flood",
                    np.arange( len(fm )),
                    [fm, ft]
                    )
    dataset.pars.scaleMonit = fm[int(low):int(high)].mean()
    dataset.pars.scaleTrans = ft[int(low):int(high)].mean()
    return dataset.pars
@}

This is the script which is run to work out the gains and offsets
of the two amplifiers. 
It works by having a scan with the shutter closed and a scan
where some region has no wafer in the beam.

@o WWMcalib.py
@{
@<imports@>
@<WWMpar@>
@<WWMdataset@>
@<WWMcalib@>

if __name__=="__main__":
    import sys
    try:
        fname = sys.argv[1]
        pars = WWMpar()
        pars.load(sys.argv[2])
        darkscan = int(sys.argv[3])
        floodscan = int(sys.argv[4])
    except:
        print "Usage: %s specfile jsonparfile darkscan floodscan"
        raise
    dataset = WWMdataset( fname, pars )
    pars = WWMcalib( dataset , darkscan, floodscan )
    print pars
    pars.save( sys.argv[2] )
@}

\subsection{ Determination of the wafer angle }

Inputs: parameters specfile listofscans

Plot signal versus angle.

Identify approximate zero (mouse click) and region to avoid due to
overflow.

Plot 1/log(signal) versus |sin(phi + phi0)|

Step 1: estimate the gradient of this to get mu t.

p = np.polyfit(abs(sin(radians(a+a0))),1/log(s),1)

OR histogram of log(signal)*|sin(phi + phi0)| -> mu t value and variance

Step 2: estimate the angular zero via:

a00 median( tan(radians(a+a0))*(1-1/mut/log(s)/abs(sin(radians(a+a0)))),",") )

since:
\[ mut/log(s) = |sin(x)| ~= |sin(x)| ( 1 + x0/tan(x) ) \]
\[ (mut/log(s)/sin(x) - 1)*tan(x) = x0 \]

Step 3: update a0 = a0 - degrees(a00)

Step 4: repeat fit removing outliers... (3 sigma on diff)

Finally, plot angle versus log(signal)*|sin(phi + phi0)| and go look for the
peaks.



\section{User interface }

This needs to be fixed. Life is short.
For now we will just hack something together using matplotlib 
and the command line.

@d UIgetrange
@{

import pylab

def UIgetrange( title, x, ylist):

    lowhigh = []
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title)
    for y in ylist:
        line, = ax.plot(x, y, ",")
    def cb(e):
        """ call back to record clicks """
        lowhigh.append(e.xdata)
        ax.plot( [e.xdata,e.xdata],ax.get_ylim(),"-")
        fig.canvas.draw()
        if len(lowhigh)==2:
            fig.close()
    cid = fig.canvas.mpl_connect('button_press_event', cb )
    pylab.show()
    return min(lowhigh[-2:]), max(lowhigh[-2:])

@}

@d UI 
@{ 
@<UIgetrange@>
@}

\section{Experimental}




\subsection{The MUSST program}

First the MUSST program, based on the one made by Ricardo Hino for testing
the rotation axis (replaced a trigger on angle step with an angular start 
position).

This generally lives in: ~blissadm/local/isg/musst/

@O wwm.mprg
@{
//************************************************************************
//************************************************************************
// MUSST program for ID11 wafer wavelength monitor
//
//  Waits for motor to reach start position
//  Records with steps constant in time
//
//   E1 = Start encoder position (moving positive)
//   DTIMER = Sampling period (beware TMRCFG defines units)
//   NPOINTS = Point number when looping over points
//   TOTPOINTS = Number of points in the scan
//               ... up to ESIZE / 4
//               ... note ESIZE is bigger when HSIZE == 0
//  
//  ALIAS are defined as it was not clear how to program without them
//
// Two programs (copy / paste ... beware!)
//  WWM_CH1:
//   ch1 motor encoder 
//   ch5 analog input
//   ch6 analog input
//  WWM_CH2:
//   ch1 motor encoder 
//   ch5 analog input
//   ch6 analog input
//************************************************************************
//************************************************************************
// point number when looping in time
UNSIGNED NPOINTS
// total numper of points in the scan
UNSIGNED TOTPOINTS = 131072
// number of encoder steps to wait before start of scan
SIGNED E1
// sampling frequency in microseconds
SIGNED DTIMER = 1000
// PROGRAM ALIASES
ALIAS A_CH1 = CH1
ALIAS A_CH2 = CH2
ALIAS A_CH5 = CH5
ALIAS A_CH6 = CH6
//**************** program starts here
PROG WWM_CH1
   // setup the timer
   CTSTOP TIMER
   CTRESET TIMER
   CTSTART ONEVENT TIMER
   // we only keep timer, encoder, ch5, ch6
   STORELIST TIMER A_CH1 A_CH5 A_CH6
   // Define the place to start the scan
   @@A_CH1 = E1
   // Wait for the event
   AT A_CH1 DO STORE
   // Now start a timescan
   NPOINTS = 1
   WHILE (NPOINTS < TOTPOINTS) DO
     @@TIMER = TIMER + DTIMER
     AT TIMER DO STORE
     NPOINTS += 1
   ENDWHILE
   EXIT NPOINTS
ENDPROG
//**************** program starts here
//**************** copy + paste of above with s/A_CH1/A_CH2/
//**************** seems easier than making it a variable
PROG WWM_CH2
   // setup the timer
   CTSTOP TIMER
   CTRESET TIMER
   CTSTART ONEVENT TIMER
   // we only keep timer, encoder, ch5, ch6
   STORELIST TIMER A_CH2 A_CH5 A_CH6
   // Define the place to start the scan
   @@A_CH2 = E1
   // Wait for the event
   AT A_CH2 DO STORE
   // Now start a timescan
   NPOINTS = 1
   WHILE (NPOINTS < TOTPOINTS) DO
     @@TIMER = TIMER + DTIMER
     AT TIMER DO STORE
     NPOINTS += 1
   ENDWHILE
   EXIT NPOINTS
ENDPROG
@}

\subsection{The SPEC scanning macros}


These macros are used to collect the data. 
They need a bit more work.

This probably should live in ~blissadm/local/spec/macros/id11

@O wwm.mac
@{

# wwm == wafer wavelength monitor

# Macro set to run a musst program for using a silicon wafer as
# wavelength monitor
#
# beam ===> diode1 =====> silicon wafer ====> diode2 
#             |                |                |
#             V                V                V
#           femto       encoded rotation      femto
#             |                |                |
#             V                V                V
#          musst_adc        musst_enc        musst_adc
#
# Record the transmission of the wafer as a function of angle
# Identify the positions of all diffractions accurately (enc. res)
# Make a fast scan
# Determine the beam wavelength and bandpass from the spectrum



# FIXME : wwm_setup macro
# FIXME : wwm_calibrate macro (eg, opening/closing main shutter)
# FIXME : globals outside of single WWM_PAR array
# FIXME : try out faster ADC card (... ordered)
# FIXME : sampling the TTL counters for V2F (RH/JMC september)
# FIXME : HSIZE set to zero on setup
# FIXME : ESIZE queried to determine memory usage (512/nchan = 131072)
# FIXME : motor channel as a parameter for the musst program
# FIXME : determination of offsets and scales semi-automatically
# FIXME : musst IO channels used to control/drive Femto amplifiers gain
# FIXME : determination of filter parameters depending on ADC card
# FIXME : is it worth doing binary data transfer for speed? ?*EDAT
# FIXME : is there any point in having stored the time?

# notes:
# /data/id11/crystallography/jon/june13/wwm_35kev/wwm_42kev_miome_speed.spc
# 
# At or below 0.4 ms speed we saw an energy vibration effect ?
# Use at least 0.8 ms ??
#
# femto gains and noise table versus time
# new musst card ordered ~ 3 kE 


global WWM_PAR[]


WWM_PAR["offset_CH5"] = 7.238e5
WWM_PAR["offset_CH6"] = -2.345e7
WWM_PAR["scale_CH5"]  = 1.0
WWM_PAR["scale_CH6"]  = 1.0
WWM_PAR["enc_chan"]  = "CH2"
WWM_PAR["mprg_name"]  = "WWM_CH2"


global WWM_MOTOR
global WWM_MOTOR_N 
WWM_MOTOR = "frz"
WWM_MOTOR_N = motor_num(WWM_MOTOR)
WWM_PAR["enc_steps_per_degree"] = motor_par(WWM_MOTOR_N,"step_size")

def wwm_set_filt(n) '{
	# Sets the filter timescale on the musst card
	# it is 2**n times the sampling (40 kHz apparently)
	# ...or could be 4 kHz, depening on the daughter card
	local two_n
	if ( n == 0 ){
		musst_comm("CHCFG CH5 ADC +-10V")
		musst_comm("CHCFG CH6 ADC +-10V")
	}
	musst_comm(sprintf("CHCFG CH5 ADC +-10V FILT %d",n))
	musst_comm(sprintf("CHCFG CH6 ADC +-10V FILT %d",n))
	p musst_comm("?CHCFG CH5")
	p musst_comm("?CHCFG CH6")
	two_n = pow(2,n)
	printf("Sampling period is 40kHz / 2^%d: %.2f Hz, %f ms",\
		n, 4.0e4/two_n, two_n/40.0)
}'

def wwm_save '{

# Saves the last collected musst scan data into the specfile

	global WWM_PAR

        local npts  i j isgname npoints nval  n 
        local long array mydata[131073][4] 
	local float array ang[131073]
	local float array signal[131073]

	nval = 4

        isgname = MUSST_AUX["default"]

 	npoints = musst_comm("?VAR NPOINTS")

	# FIXME: binary transfer? Save as edf?
        n = musst_getdata(isgname, npoints, nval, mydata)

        # This is from _head in ascan
        #
        if (DATAFILE != "") {
                local i,j,z,s
                ond; offt

               HEADING = sprintf("wwscan %g %g %d %g",\
			 WWM_PAR["start"],\
			 WWM_PAR["end"],\
			 WWM_PAR["nstep"],\
			 WWM_PAR["ctime"])

                printf("\n#S %d  %s\n#D %s\n",++SCAN_N,HEADING,DATE)

                if (_ctime < 0) printf("#M %g  (%s)\n", -_ctime,cnt_name(MON))
                else printf("#T %g  (%s)\n", _ctime, cnt_name(sec))
                _head_par G 0
                _head_par U 1
                _head_par UB 3
                _head_par Q 4
                printf("#Q %s\n", _hkl_val)
                for (i=0; i<MOTORS; i+= 8) {
                        s = sprintf("#P%d ",i/8)
                        for (j=i; j<i+8 && j<MOTORS;) {
                                if (motor_name(mA[j]) != "unused")
                                        s = s sprintf("%.8g", A[mA[j]])
                                if (j%8 == 7)
                                        break
                                s = s " "
                                j++
                        }
                        print s

                }
                Fheader
                user_Fheader
                z = _ctime < 0? sec:MON

        # MUSST specfic code starts here!
        # 
        # number of counters
        printf("#N %d\n", 6 )
        # names of counters
	#            0      1    2    3
        printf("#L  angle  Timer  %s  CH5  CH6  signal\n",WWM_PAR["enc_chan"])

	# And now the data
        for(i=0; i< npoints ; i++){
		signal[i] =  mydata[i][2]/mydata[i][3]
		ang[i] =  mydata[i][1]/WWM_PAR["enc_steps_per_degree"]
		printf("%f  ",ang[i])
                for(j=0; j<nval ; j++){
                        printf("%d  ",mydata[i][j])
                        }
		printf("%f  ",signal[i])

                printf("\n")
                }
        offd ;  ont
        }

	wwm_plot( ang[:npoints-1], signal[:npoints-1] )

}'


def wwm_plot( x, y )'{
	plot_cntl("erase")
	plot_range("auto","auto","auto","auto")
	array_plot(x, y )
}'

def wwm_move(_pos)'{
 global WWM_MOTOR_N
 waitmove
 getangles
 A[WWM_MOTOR_N]=_pos
 move_em
 _update(WWM_MOTOR)
}'

def wwm_scan '{

#FIXME!!!


}'

def _wwm_scan(start, end, nstep, ctime)'{ 
 global WWM_MOTOR WWM_MOTOR_N
 local totaltime totalsteps velocity stime nfilt 
 local start_pos_steps undershoot

 if(end <= (start + 0.1)){
   printf("You have to scan in a positive sense\n")
   exit
 }
 
 if( nstep > 131072 ){
   printf("Too many steps\n")
   exit
 }

 WWM_PAR["start"] = start
 WWM_PAR["end"] = end
 WWM_PAR["nstep"] = nstep
 WWM_PAR["ctime"] = ctime


 totaltime = nstep * ctime
 totalsteps = (end - start)*motor_par(WWM_MOTOR_N, "step_size")

 # compute motor speed required 	
 velocity = totalsteps/totaltime 

 if( velocity > motor_par( WWM_MOTOR_N, "config_velocity" )){
   printf("Your scan is too fast, try again\n")
   exit
 }

 printf("Moving %s to start position %f\n",WWM_MOTOR,start)

 # Page of integral calculations shows that the undershoot
 # comes out as the average speed while accelerating times 
 # the acceleration time 
 
 # convert from milliseconds
 acctime =  motor_par(WWM_MOTOR_N, "acceleration")/1e3 
 # average speed
 avgaccvel = (motor_par(WWM_MOTOR_N, "base_rate") + velocity)*0.5
 # steps to reach constant speed
 undershoot = avgaccvel * acctime / motor_par(WWM_MOTOR_N, "step_size")
 printf("Undershoot should be %f",undershoot)

 wwm_move( start - undershoot )
 # set speed
 printf("Totaltime= %f, totalsteps= %f, velocity= %f\n",\
	totaltime,totalsteps,velocity)
 motor_par( WWM_MOTOR, "velocity", velocity)
 printf("Set motor speed to %f\n", motor_par( WWM_MOTOR_N, "velocity"))
 
 # Set up and run the musst program
 musst_comm(sprintf("VAR TOTPOINTS %d",nstep))
 
 # FIXME: the 1e6 comes from TMRCFG==1MHZ in musst
 musst_comm(sprintf("VAR DTIMER %d",ctime*1e6))

 
 # Set the start angle position
 musst_comm(sprintf("VAR E1 %d",start*motor_par(WWM_MOTOR_N, "step_size")))

 # Set the filter time to smooth over the ctime
 nfilt =  int(log(40000*ctime)/log(2))
 wwm_set_filt( nfilt )

 stime = time()

 # run the musst program
 musst_comm(sprintf("RUN %s",WWM_PAR["mprg_name"]))

 printf("Scanning the wafer now!\n")
 # Move the motor to the end position 
 wwm_move(end+undershoot)

 printf("Motor move finished in %.2f /s\n",time()-stime)
 # set speed back to normal
 motor_par( WWM_MOTOR, "velocity", motor_par(WWM_MOTOR,"config_velocity"))

 while( musst_comm("?STATE") != "IDLE" ){
	printf("Waiting for musst %f",time()-stime)
	sleep(1)
 }

 wwm_save
 printf("All done, back to you %.2f /s\n",time()-stime)

}'

def _wwm_collect_full(ctime) '{

local st
# scan near 90 to get transmission
exit
shclose
# Offsets
_wwm_scan(-100,-79.8,60200,ctime)
shopen
# Linearity - beam will shoot past
smvice hx 6.68
_wwm_scan(-100,-79.8,60200,ctime)
# Now scan with crystal on the axis
smvice hx 1.68
for(st=-180;st<161;st+=20){
 _wwm_scan(st,st+20.2,60200,ctime)
 }

}'
 
def _wwm_collect_mini(ctime) '{

local st
_wwm_scan(20,42.5,131070,ctime)
}'

def _wwm_frz_full(ctime) '{

for(st=0;st<359;st+=22.5){
  _wwm_scan(st,st+22.5,131070,ctime)
 }

}'

def wwm_monitor_drift '{

while(1){
remote_async("lid112:optics","lrock 0.01")
wait(0x08)
fprintf(DATAFILE,"#C rz1 %f\n",remote_eval("lid112:optics","A[rz1]"))
fprintf(DATAFILE,"#C rz2 %f\n",remote_eval("lid112:optics","A[rz2]"))
_wwm_collect_mini(0.0016)
#sleep(60*10)
}

}'
@}



\section{Building the code}

@o bldwwm.bat
@{
nuweb wwm
pdflatex wwm
@}

\section{ Literature }

The Indexing of Single-Crystal X-ray Rotation Photographs.
J. R. Milch and T. C. Minor.
J. Appl. Cryst. (1974) 7, 502

Fable geometry document...

Bond (196?)

Determination of X-ray flux using silicon pin diodes
R. L. Owen, J. M. Holton, C. Shulze-Briese and E. F. Garman.
J. Synch. Rad. (2008) 16, 143-151.

Dynamical Diffraction of X Rays by Perfect Crystals
Rev Mod Phys. 36(3) 681-717 (1964)
B. W. Batterman and H. Cole

\section{ Smoothing code from the net }

This comes from the scipy wiki (http://wiki.scipy.org/Cookbook/SavitzkyGolay)

@d savitzky
@{
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
@}

\section{ Unsorted codes }

@d unsorted
@{
def find_peaks( ang, sig, npxmin = 10 ):
    """
    Find connected groups of pixels...
    """
    from ImageD11.connectedpixels import blobproperties, connectedpixels,\
        blob_moments, s_1, s_I, mx_I
    import numpy as np
    # median is the background
    sm = np.median(sig)
    # sigma level
    st = sig.std()
    # norman is our normalised array, 2D, nx1 for peaksearching
    norman = ((sig - sm)/st).astype(np.float32)
    anorman = (ang*norman).astype(np.float32)
    norman.shape = norman.shape[0],1
    anorman.shape = anorman.shape[0],1
    # label array
    b = np.zeros( norman.shape, np.int)
    print b.shape
    npks = connectedpixels( norman, b, 1.0, 0 )
    print npks
    #
    sig_pks = blobproperties( norman, b, npks, 0)
    # we want:
    #  width
    allwidth = sig_pks[:,s_1]
    # Filter out anything with less than npxmin pixels
    msk = allwidth >= npxmin
    widths = np.compress( msk, allwidth)
    #  centroid in angle (sum ints*a / sum(ints)
    print npks
    ang_pks = blobproperties( anorman , b, npks, 0)
    centroids = np.compress( msk, ang_pks[:,s_I]/sig_pks[:,s_I])
    #  area or height ?
    areas = np.compress( msk, sig_pks[:,s_I]*st )
    heights = np.compress( msk, sig_pks[:,mx_I]*st )
    return centroids, areas, widths, heights
    

if __name__=="__main__":
    data = readspc("wwm_42kev/wwm_42kev_0.0016.spc")
    ang,sig = makesig(data)
    c,a,w,h = find_peaks( ang, sig)

    order = np.argsort(-h)

    c[order[:8]]
#Out[150]: 
#array([ 147.47907659,   38.00419486,  -32.52035246, -141.99559782,
#       -147.3868076 ,  142.08722053,   32.6130885 ,  -37.91206404])
#UBI:
#5.29317363 -0.0905953467 -0.0683694637
#-0.00740113112 3.77238274 -3.81213226
#0.100616448 3.71563154 3.88134886


@}


\end{document}
