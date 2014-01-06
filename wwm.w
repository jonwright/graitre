% -*- indent-tabs-mode:nil -*-
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

\section{TODO}

- Number of peaks versus exposure time is not constant when using
 the sigma threshold. More peaks for shorter time.
- improve mut fitting (single lsq seems best)
- single gui to open file, calib offsets, fit mut and peaksearch.
- peak assignment to hkls




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
# standard library
import math, json, pprint, time

# standard science stack
import numpy as np, pylab

# Fable code from Jon
from ImageD11 import unitcell, gv_general, transform, grain
from ImageD11.connectedpixels import blobproperties, connectedpixels,\
      blob_moments, s_1, s_I, mx_I, bb_mn_f, bb_mx_f
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

def tandeg(x):
    """ tangent of angle in degrees """
    return np.tan(np.radians(x))

@<rotationmatrices@>

@<atan2dvecs@>  


@}

A quick bit of code for finding the angle between a pair of vectors in 2D,
and also the derivatives. 
Just to explain the cross product part in there. 
We want the angle between the two vectors $k_{xy}$ and $g_{xy}$
which are in the $xy$ plane. 
We consider them as 3D vectors $k_x,k_y,0$ and $g_x,g_y,0$.
Now the cross product is up the z axis and the length
is $k_x.g_y - k_y.g_x$ which is equal to $|k||g|sin(\theta)$.
We then mix this with the sin and atan2 to get the angle we
want
Due to laziness we call upon sympy to check the derivatives.

@o vec2anglesdiff.py
@{
from sympy import Symbol, diff, atan2, simplify
a,b,c,d = [Symbol(x) for x in 'abcd']
e = atan2( a*d - b*c, a*c + b*d) 
print e
for s in a,b,c,d:
    print s,":",simplify(diff(e,s))
@}

The output of the script is:

\begin{verbatim}
atan2(a*d - b*c, a*c + b*d)
a : b/(a**2 + b**2)
b : -a/(a**2 + b**2)
c : -d/(c**2 + d**2)
d : c/(c**2 + d**2)
\end{verbatim}

We use this below:

@d atan2dvecs 
@{
def anglevecs2D(a, b, c, d):
    return np.arctan2( a*d - b*c, a*c + b*d) 
    
def derivanglevecs2d( a, b, c, d, ab2=None, cd2=None):
    if ab is None: 
        ab2 = a*a+b*b
    if cd is None:
        cd2 = c*c+d*d
    return {
            'a':  b / ab,
            'b': -a / ab,
            'c': -d / cd,
            'd':  c / cd,
            }
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

def cleanbuf( buf, nvals ):
    for line in buf:
        if line[0] == "#":
            continue
        vals = line.split()
        if len(vals) == nvals:
            yield [float(v) for v in vals]

def decodescan( buf, ncols ):
    return np.array( [row for row in cleanbuf(buf, ncols)])


def readspc(fname):
    """ minimalist spec file reading """
    import cStringIO
    scans = []
    coltitles = []
    start = time.time()
    with open(fname) as f:
        lines = f.readlines()
        scanlines = []
        for i,line in enumerate(lines):
            if line[0:2] == "#S":
                scanlines.append(i)
            if line[0:2] == "#L":
                titles = line[2:].strip().split("  ") 
                coltitles.append(titles)
        ns = len(scanlines)
        for i in range(ns):
            s = scanlines[i]
            if i+1 < ns:
                e = scanlines[i+1]
            else:
                e = len(lines)
            scans.append( lines[s:e] ) 
        print "time to chop scans into buffers",fname,time.time()-start
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
    'zeroAngle'  : -90.0,
    'scaleMonit' : 1.0,
    'scaleTrans' : 1.0,
    'scaleAngle' : 160000.0,
    'chanMonit'  : 'CH6',
    'chanTrans'  : 'CH5',
    'chanAngle'    : 'CH1',
    'mut'   : 0.05,
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

@O wwmfrz.json
@{
{
    "chanAngle": "CH2",
    "chanMonit": "CH5",
    "chanTrans": "CH6",
    "mut": 0.017244892180296886,
    "scaleAngle": 17493.333330000001,
    "scaleMonit": -189161910.72070956,
    "scaleTrans": -167473298.25006324,
    "zeroAngle": -0.058586303894149205,
    "zeroMonit": 655356.36976666667,
    "zeroTrans": 471232.09580000001
}
@}

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

\subsection{ Refactor that to make a crystal }

We recycle the WWMpar code via inheritance to make a crystal object:

@d WWMcrystal
@{
class WWMcrystal( WWMpar ):
    a = 5.43094 # silicon, value to be debated.
    dmin = None
    def set_wvln(self, wvln):
        self.wvln = wvln
        self.bmatrix = [ [ wvln/self.a, 0, 0],
                         [ 0,  wvln/self.a, 0],
                         [ 0,  0, wvln/self.a] ]

    def generate_hkls( self ):
        """
        Generate all of the hkls which can be reached 
        a is the cubic unit cell parameter
        """
        uc = unitcell.unitcell( [ self.a,self.a,self.a,90,90,90],"F" )
        if self.dmin is None:
            uc.gethkls_xfab( 2.0/self.wvln, "Fd-3m" ) 
        else:
            uc.gethkls_xfab( 1.0/self.dmin, "Fd-3m" ) 
        self.hkls = np.array([p[1] for p in uc.peaks])
        self.readfhkl( self.fhklfile )
        self.icalc = np.array( [self.getfhkl(*h) for h in self.hkls] )
        return self.hkls

    def set_orientation( self, umatrix):
        """ 3x3 list please """
        self.umatrix = umatrix

@< generategve @>
@< crystalderivatives @>
@< crystalintensities @>

    def calcall( self ):
        self.generate_hkls()
        self.generate_gve()
        self.generate_derivatives()


def maketestcrystal():
    c = WWMcrystal()
    c.set_wvln( 1.0 )
    u0= np.dot( np.dot(  rotmatx(0.1),rotmaty(0.2)),rotmatz(0.3))
    c.set_orientation( u0 )
    c.calcall()
    return c

@}

Now we get to the fun parts, computing the peak positions and their
derivatives. 
We make this the job of the crystal object.

First we will compute the g-vectors corresponding to the hkls and u (and the
b matrix). 
The bmatrix here carries the wavelength term as well as cell parameter.

@d generategve
@{
    def generate_gve(self):
        assert len(self.hkls)>1
        self.ub = np.dot( self.umatrix, self.bmatrix )
        self.gve = np.dot( self.ub, self.hkls.T )
@}
        
And now some derivatives. 
These are for the gvectors. 
Surprisingly painless.

@d crystalderivatives 
@{
    VARIABLES = ['wvln', 'rotx', 'roty', 'rotz']
    def generate_derivatives(self):
        self.derivatives = {}
        # wavelength 
        # db / dwvln
        self.derivatives['wvln'] = self.gve / self.wvln
        self.derivatives['rotx'] = np.cross( self.gve.T, [1,0,0] ).T
        self.derivatives['roty'] = np.cross( self.gve.T, [0,1,0] ).T
        self.derivatives['rotz'] = np.cross( self.gve.T, [0,0,1] ).T
@}

Functions to make the rotation matrices which correspond
to rotations around x, y and z

@d rotationmatrices
@{
def rotmatx(r):
    """ r is in radians """
    c = np.cos(r)
    s = np.sin(r)
    return [ [ 1,0,0], [0, c, s], [0,-s,c] ]
def rotmaty(r):
    """ r is in radians """
    c = np.cos(r)
    s = np.sin(r)
    return [ [ c,0,-s], [0, 1, 0], [s,0,c] ]
def rotmatz(r):
    """ r is in radians """
    c = np.cos(r)
    s = np.sin(r)
    return [ [ c,s,0], [-s, c, 0], [0,0,1] ]

def getrxryrz(u):
    """ 
    Converts a U matrix to give rx, ry, rz angles
    Wikipedia page X1Y2Z3 with Tait-Bryan angles \
    """
    c2s1 = -u[1,2]
    c2c1 = u[2,2]
    r1   = -np.arctan2(c2s1,c2c1)
    c2c3 = u[0,0]
    c2s3 = -u[0,1]
    r3 = -np.arctan2( c2s3, c2c3 )
    s2 = u[0,2]
    if abs(np.sin(r3)) > 0.5:
        c2 = c2s3 / np.sin(r3)
    else:
        c2 = c2c3 / np.cos(r3)
    r2 = -np.arctan2( s2, c2 )
    if 1:
        utest = np.dot(np.dot(rotmatx(r1),rotmaty(r2)),rotmatz(r3))
        assert abs(utest-u).ravel().sum() < 1e-10
    return r1,r2,r3

@}

\subsection{Computed intensities}

To get the computed F(hkl) values for peaks is quite a lot of 
work. 
A quick workaround is to read the output from Lazy-Pulverix
on the ICSD web server to get values.

@d crystalintensities
@{
    fhklfile = "si.fhkl" # from icsd/lazy
    def readfhkl(self,fname):
        """ reads icsd/lazy file """
        self.fhkl = {}
        for line in open("si.fhkl").readlines():
            items = line.split()
            if len(items)<5:
                continue
            try:
                hkl = [int(v) for v in items[:3]]
            except:
                continue
            intensity = float(items[13])**2 # A**2 ok for silicon only
            self.fhkl[tuple(hkl)] = intensity

    def getfhkl(self,h,k,l):
        """ returns I(hkl) = F(hkl)*F'(hkl) given hkl """
        pos = [abs(h),abs(k),abs(l)]
        pos.sort(reverse=True)
        pos = tuple(pos)
        if self.fhkl.has_key(pos):
            return self.fhkl[pos]
        else:
            #        print "miss",h,k,l,pos
            return 0.0

    def filter_hkls( self, fracimax):
        """
        remove hkls which are below some fraction of the strongest
        """
        imax = self.icalc.max()
        mask = self.icalc/imax > fracimax
        print self.hkls.shape
        self.hkls = np.compress( mask, self.hkls, axis=0 )
        print self.hkls.shape
        self.icalc = np.compress( mask, self.icalc )
@}


\subsection{Crystal testcases}

Finally a testcase for the WWMcrystal, to check we have the rotation
matrix definitions the right way around too.

@o testcrystal.py
@{
@<imports@>
@<utilities@>
@<WWMpar@>
@<WWMcrystal@>




def testwvln():
    print "test wvln"
    c = maketestcrystal()
    # now test the derivatives are OK or not
    g0 = c.gve.copy()
    c.set_wvln( 1.001 )
    c.generate_gve()
    dgdw = (c.gve - g0)/0.001
    err = abs( dgdw - c.derivatives['wvln'] ).ravel()
    iworst = err.argmax()
    assert err[iworst] < 1e-6
    #  print err[iworst],dgdw[:,iworst/3]-c.derivatives['wvln'][:,iworst/3],
    # print c.derivatives['wvln'][:,iworst/3]

def testrot(axis):
    print "testrot",axis
    c = maketestcrystal()
    # now test the derivatives are OK or not
    g0 = c.gve.copy()
    # rotation around the x axis
    dx = 1e-6
    if axis == 'x': dr = rotmatx(dx)
    if axis == 'y': dr = rotmaty(dx)
    if axis == 'z': dr = rotmatz(dx)
    c.set_orientation( np.dot(dr, c.umatrix ) )
    c.generate_gve()
    dgdr = (c.gve - g0)/dx
    key = 'rot'+axis
    err = abs( dgdr - c.derivatives[key] )
#    for i in range(0,910,1):
#    	print i,dgdr[:,i],err[:,i]
    iworst = err.ravel().argmax()
#    print err.ravel()[iworst]
    assert err.ravel()[iworst] < 1e-6
if __name__=="__main__":
   testwvln()
   for axis in ('x','y','z'):
      testrot(axis)

@}


\subsection{ Diffractometer class }

This takes a crystal (list of hkls, perhaps intensity) and uses an orientation
to get computed omega angles and derivatives.

@d WWMdiffractometer
@{
class WWMdiffractometer( object ):
    def __init__(self, crystal, pars, axistilt=0):
        """
        crystal is a WWMcrystal instance (gives hkl list)
        u is a 3x3 orthogonal rotation matrix (z on axis)
        axistilt sin(radians) is the sin of the angle between the axis and
                 the plane perpendicular to the beam
        pars : WWMparameters json file
        """
        self.crystal = crystal
        self.axistilt = axistilt
        # For crystal normal
        self.pars = WWMpar( pars )

@<computek@>
@<computeomegas@>
@<pathlengths@>

    def calcall(self):
        self.computek()
        self.computeomegas()

@<derivDiffractometer@>

@}   

The first part is to compute the k vectors for a given crystal (eg, the 
scattering vectors at the moment diffraction occurs).
Not all peaks can be measured, some will be along the axis. 
These are noted by self.mask.
From the Milch and Minor magic:

@d computek
@{
    def computek(self):
        """ This is from Milch and Minor 1974 """
        g = self.crystal.gve
        h = self.crystal.hkls
        # Make two copies
        self.g = np.concatenate((g,g),axis=1)
        self.h = np.concatenate((h,h),axis=0)
        #  print "hkls shape",self.h.shape,"gve.shape",self.g.shape
        self.fri = np.concatenate((np.ones(g.shape[1]),-np.ones(g.shape[1])))
        self.modg2 = (self.g*self.g).sum(axis=0)
        self.kz = self.g[2,:]
        num = - (self.modg2 - 2*self.axistilt*self.kz)
        den = 2*np.sqrt( 1 - self.axistilt*self.axistilt )
        # Normally axistilt should be a smallish number
        self.kx = num / den
        # Some peaks never diffract as they are up the axis
        arg = self.modg2 - self.kx*self.kx - self.kz*self.kz
        self.mask = arg < 0
        # Allow sqrt negative as nan for now
        # self.fri is the sign to take on the Friedel pair
        self.ky = np.sqrt( arg ) * self.fri
        self.icalc = np.concatenate( (self.crystal.icalc, self.crystal.icalc ) )
@}


And now for the angles, given the k and g vectors:

@d computeomegas
@{
    def computeomegas(self):
        g = self.g
        gx = g[0,:]
        gy = g[1,:]
        gz = g[2,:]
        self.omega  = anglevecs2D( gx, gy, self.kx, self.ky)

@}

To estimate the peak intensities we will need to distinguish 
between the case of Laue diffraction (transmission through the 
crystal) and Bragg reflection (diffraction comes out of the 
same side of the crystal.
For a given omega angle we will need to know the direction of
the surface normal of the crystal.
We should get this from the mut computation.

We have the beam in as $k_{in}$, the beam out as $k_{out}$ 
and the scattering vector as $k = k_{out} - k_{in}$

@d pathlengths
@{
    def pathlengths(self):
        """
        The ratio of the output pathlength for the diffracted ray
        in comparison to the input pathlength
        At omega == self.pars.zeroAngle the wafer is parallel to the beam
        """
        a0 = self.pars.zeroAngle # degrees
        # along the beam, sin zero is zero
        n = np.array( (np.sin( self.omega) ,
                       np.cos( self.omega) ,
                       np.zeros( len(self.omega))  )     )
        print "nshape",n.shape
        # Path on input direction is roughly nx (beam on 1,0,0)
        # cosine of angle between normal and (1,0,0)
        pin = 1.0/n[0,:]
        # zn is zero, we assume it is roughly on axis here (not quite true)
        # FIXME - intensity refinement will need the fit the normal direction
        #
        # Output direction
        ko = np.array((self.kx + 1 , self.ky, self.kz)) 
        normedko = ko / np.sqrt((ko*ko).sum(axis=0))
        print "koshape",ko.shape
        #
        # cosine with output direction
        cospko = (normedko*n).sum(axis=0) 
        po = 1.0/cospko
        print "poshape",po.shape
        #
        # Both of these are signed quantities. 
        # Same sign == Bragg
        # opposite sign == Laue
        #
        print "pinshape",pin.shape
        for i in range(10):
            print self.omega[i], pin[i], po[i],
            print self.kx[i],self.ky[i],self.kz[i]
        #

        # plr == path length ratio
        self.plr = np.arctan2(pin, po)
        print self.plr.shape
        return self.plr
@}

        
And a quick testcase for the diffractometer stuff:

@d testdiff1
@{
def testdiff1():
    d = WWMdiffractometer( maketestcrystal(), testpars )
    d.calcall()
    # print d.h.shape
    for i in range(50):
        print "%d "%(i),"%4d %4d %4d"%tuple(d.h[i]),
        print "%.5f %.5f %.5f "%tuple(d.g[:,i]),
        print "%.5f "%(d.modg2[i]),
        print "%.5f %.5f %.5f "%(d.kx[i],d.ky[i],d.kz[i]),
        print "%.5f "%(d.omega[i]),
        print d.mask[i]
@}


@o testdiff.py
@{

@<imports@>
@<utilities@>
@<WWMpar@>
@<WWMcrystal@>
@<WWMdiffractometer@>
@<testdiff1@>
if __name__=="__main__":
   import sys
   testpars=sys.argv[1]
   testdiff1()
@}

\subsection{Angle derivatives }

Finally, the part we have all been waiting for.
Derivatives of the angles of the hkl reflections with respect to the 
variables.


@d derivDiffractometer
@{
    # ['wvln', 'rotx', 'roty', 'rotz']
    VARIABLES=  WWMcrystal.VARIABLES + [
       'axistilt' ]
    def generate_derivatives(self):
        # g = self.crystal.gve
        # self.modg2 = (g*g).sum(axis=0)
        # we have dg/dv in self.crystal.derivatives
        dg = self.crystal.derivatives
        # dict comprehensions, python 2.7+
        dkzdv = dict( (v,  dg[v][:,2]) 
                for v  in WWMcrystal.VARIABLES )
        dmodg2dv = dict( (v, (2 * g * dg[v]).sum(axis=0) )
                for v in WWMcrystal.VARIABLES  )
        # num = - (self.modg2 - 2*self.axistilt*self.kz)
        dnumdv = - ( dmodg2dv ) * FIXMEHERE
        # den = 2*np.sqrt( 1 - self.axistilt*self.axistilt )
        raise "not implemented yet"
        


               

@}


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

@d WWM
@{
@<imports@>
@<utilities@>
@<UI@>
@<readspc@>

@<WWMpar@>

class WWMdataset(object):
    def __init__(self, fname, pars):
        self.fname = fname
        self.scans, self.coltitles = readspc( fname ) 
        self.pars = pars
    def checkscan(self, scannumber):
        """ check scan is indeed in the file
            and decodes the ascii on demand """
        assert scannumber >= 0 and scannumber < len(self.scans), \
                    (s,scannumber, len(self.scans))
        if not isinstance(self.scans[scannumber], np.ndarray):
            self.scans[scannumber] = decodescan(self.scans[scannumber], 
                                        len(self.coltitles[scannumber]))
    def getcol(self, scan, name):
        """ Accepts a single scan or a list of scans """
        if np.iterable(scan): 
            # Bounds check first
            i = self.coltitles[scan[0]].index(name)
            for s in scan:
                self.checkscan(s)
            assert i == self.coltitles[s].index(name)
            return np.concatenate( [self.scans[s][:,i] for s in scan] )
        else:
            if name in self.coltitles[scan]:
                i = self.coltitles[scan].index(name)
            else:
                print name
                print self.coltitles[scan]
                raise
            self.checkscan(scan)
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
# we make the calib a class method for a dataset
@<WWMcalib@>
# And we have a class method for fitting the mut and wafer zero:
@<WWMmut@>
# And we have a class method for computing the profile
@<WWMabsprof@>
@<WWMpeaksearch@>


@}

\subsection{Dark and Gain Calibration}

We combine a set of experimental parameters with a scan to 
get to the calibration. 
Minimal inputs are the specfilename and scan numbers, we can also 
add the scan motor steps and angle channel.

@d WWMcalib
@{

    def calib( dataset, darkscan, floodscan ):
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
@<WWM@>
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
    pars = dataset.calib( darkscan, floodscan )
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

@d WWMmut
@{
    def mut(self, scans, cutrange=2.5, debug=True):
        """
        a0 = zero angle for wafer, ask if None
        cutrange = excluded angular range
        """
        a = self.getAngle(scans)
        s = self.getSignal(scans)
        a0 = self.pars.zeroAngle
        if a0 is "None":
            a0,junk = UIgetxy( "Zero angle for wafer please", a,[s,])
            print "Got zero angle as",a0,junk
        acor = angmoddeg(a-a0)
        msk = (abs(acor) > cutrange)&(abs(acor)<180-cutrange)
        af = np.compress( msk, a )
        sf = 1.0/np.log(np.compress( msk, s )) # log signal here
        ncycle = 5
        for i in range(ncycle):
            x = abssindeg(af-a0) # changes per run
            p = np.polyfit( x, sf, 1 ) # offset term ??
            mut = -1/p[0]
            a0a = tandeg(af-a0)*(sf*mut/x + 1)
            a00 = np.median(a0a)
            if debug or i == ncycle-1: 
                print "mut is",mut, p
                pylab.figure()
                pylab.plot(x*np.sign(af-a0), sf,",")
                print "a00 is",a00, np.degrees(a00),a0
                print "deviations",a0a.mean(),a0a.std()
                pylab.figure()
                pylab.plot( af-a0, a0a ,",")
                pylab.show()
                raw_input()

            a0 += np.degrees(a00)
            if i < ncycle-1:
                outliers = abs(a0a - a00)/a0a.std() < 3
                print "Suppressed",len(outliers)-outliers.sum()
                af = np.compress( outliers, af )
                sf = np.compress( outliers, sf )
        self.pars.mut = mut
        self.pars.zeroAngle = a0
@}      

@o WWMmut.py
@{
@<WWM@>
if __name__=="__main__":
    import sys
    try:
        fname = sys.argv[1]
        pars = WWMpar()
        pars.load(sys.argv[2])
        scans = [int(v) for v in sys.argv[3:]]
    except:
        print "Usage: %s specfile jsonparfile scans"
        raise
    dataset = WWMdataset( fname, pars )
    dataset.mut( scans, debug = False )
    dataset.pars.save(sys.argv[2])
@}

\subsection{Peak searching}

Assuming we have all calibrated we want to correct the data 
for the mut and go look for peaks.

@d WWMabsprof
@{
    def absprof(self, scans):
        """
        absorbtion profile.
        Takes log(signal) / abs(sin(angle))
        """
        a = self.getAngle(scans)
        s = self.getSignal(scans)
        a0 = self.pars.zeroAngle
        self.angle = a-a0
        x = abssindeg(self.angle) # changes per run
        mut = self.pars.mut
        corr = -mut-np.log(s)*x
        self.corr = corr
@}      


@d WWMpeaksearch
@{
    def peaksearch(self, scans, filename, npxmin=10, plot=True):
        """ 
        finds peaks in self.corr

        """
        if not hasattr(self,"corr"):
            self.absprof( scans )
        # median is the background
        ang = self.angle
        sm = np.median(self.corr)
        # sigma level
        st = self.corr.std()
        threshold = sm+st # 1 sigma
        pcorr = self.corr
        acorr = ang*self.corr
        print self.corr.shape, ang.shape
        a2corr = ang*acorr
        pcorr.shape = pcorr.shape[0],1
        acorr.shape = pcorr.shape[0],1
        a2corr.shape = pcorr.shape[0],1
        # label array
        b = np.zeros( pcorr.shape, np.int)
        print "Searching in ",b.shape,
        npks = connectedpixels( pcorr, b, threshold, 0 )
        print "found",npks,"peaks"
        #
        sig_pks = blobproperties( pcorr, b, npks, 0)
        ang_pks = blobproperties( acorr , b, npks, 0)
        ang2_pks = blobproperties( a2corr , b, npks, 0)
        # we want:
        #  npx
        npx = sig_pks[:,s_1]
        # Filter out anything with less than npxmin pixels
        msk = npx >= npxmin
        #  centroid in angle (sum ints*a / sum(ints)
        print "After removing zingers",npks
        centroids = np.compress( msk, ang_pks[:,s_I]/sig_pks[:,s_I])
        #  area or height ?
        areas = np.compress( msk, sig_pks[:,s_I] )
        heights = np.compress( msk, sig_pks[:,mx_I] )
        # width 
        low = np.compress( msk, sig_pks[:,bb_mn_f]).astype(int)
        high = np.compress( msk, sig_pks[:,bb_mx_f]).astype(int)
        widths = ang[high] - ang[low]
        # sigma
        # sum( x*x*I )
        xxI = np.compress( msk, ang2_pks[:,s_I])
        vari = xxI / areas - centroids*centroids
        stof = 8*np.log(2)
        sigma = np.sqrt( abs( vari )*stof )
        if plot:
            pylab.plot(self.angle,self.corr,",")
            pylab.plot(centroids, heights,"+")
            pylab.plot(centroids-widths/2, heights*0.3,"|")   
            pylab.plot(centroids+widths/2, heights*0.3,"|")   
            pylab.plot(centroids-sigma/2, heights*0.5,"|")   
            pylab.plot(centroids+sigma/2, heights*0.5,"|")   
            pylab.show()
        self.pars.centroids = list(centroids)
        self.pars.areas = list(areas)
        self.pars.widths = list(widths)
        self.pars.heights = list(heights)
        self.pars.save( filename )    
        return centroids, areas, widths, heights
@}


@o WWMpeaksearch.py
@{
@<WWM@>
if __name__=="__main__":
    import sys
    try:
        fname = sys.argv[1]
        pars = WWMpar()
        pars.load(sys.argv[2])
        outfile = sys.argv[3]
        scans = [int(v) for v in sys.argv[4:]]
    except:
        print "Usage: %s specfile jsonparfile scans"
        raise
    dataset = WWMdataset( fname, pars )
    dataset.peaksearch(scans, outfile)
@}


\section{Assignment of peaks}

Given two lists of peaks try to match up which ones should go together.
This can be observed and computed peaks or also peaks from two 
different samples.
We will assume both sets of peaks come in sorted order.

@d matchpeaks
@{

def matchpeaks(x1, x2, tol=0.1):
    """
    Match up two peaks lists x1 and x1 assuming both are sorted
    all matches within tol are returned
    """
    i1 = 0
    i2 = 0
    matches = []
    while i1<len(x1) and i2<len(x2):
        v1 = x1[i1]
        v2 = x2[i2] 
        if abs(v1-v2) < tol:
            matches.append( (i1, i2) )
        if v1>v2: 
            i2 += 1
        else:
            i1 += 1        
    return matches
@}

@O WWMpeakassign.py
@{
@<imports@>
@<utilities@>
@<WWM@>
@<WWMpar@>
@<WWMcrystal@>
@<WWMdiffractometer@>


from ImageD11.simplex import Simplex

class simplexfitter(object):
    def __init__(self,energy,ubifile,peaksfile):
        c = WWMcrystal()
        c.set_wvln(12.3985/energy)
        # c.dmin = c.a/np.sqrt(3)*0.9
        gr = grain.read_grain_file( ubifile )[0]
        u0 = gr.u
        rx, ry, rz = getrxryrz( u0 )
        self.start = [rx,ry,rz]
        utest= np.dot( np.dot(  rotmatx(rx),rotmaty(ry)),rotmatz(rz))
        if not abs(utest -u0).sum() < 1e-10:
           print u0, utest
           raise
        self.start = [rx, ry, rz]
        c.set_orientation(u0)
        c.generate_hkls()
        c.filter_hkls( 0.01 )
        c.generate_gve(  )
        d = WWMdiffractometer(c, peaksfile)

        d.computek()
        print len(d.icalc)
        d.computeomegas()
        p = WWMpar(peaksfile)
        p.heights = np.array(p.heights)
        obs = np.array(p.centroids)+p.zeroAngle
        self.obspeaks=obs
        # for each calc peak find the closest obs peak
        pairs = []
        oc = np.degrees(d.omega)
        for obs_j in range(len(obs)):
            diffs = abs( obs[obs_j] - oc  )
            inds = np.argsort( diffs )
            i = inds[0]
            #print diffs[i]
            #1/0
            #     print oc, obs_j, diffs[obs_j]
            if diffs[i] < 0.15 and diffs[inds[1]]>2*diffs[inds[0]]:
               #print "got"
               pairs.append((i,obs_j))
        print len(pairs)
        oms =  np.array( [ (np.degrees(d.omega[i]), obs[j])        for i,j in pairs ])
        print oms.shape
        ints = np.array( [ (d.icalc[i], p.heights[j])        for i,j in pairs])
        pathlenrat = np.take( d.pathlengths(), [i for i,j in pairs])
        if 0:
            pylab.figure()
            pylab.plot( oms[:,0], oms[:,0]-oms[:,1],"o")
            pylab.figure()
            pylab.plot( oms[:,0], ints[:,0]/ints[:,1],"o")
            pylab.figure()
            pylab.plot( pathlenrat, ints[:,0]/ints[:,1],"o")
            pylab.figure()

        dataset = WWMdataset( "crystal_ori.spc", p )
        dataset.peaksearch([7,8,9,10,11,12], "crystal_ori.spc.json")
        dom = np.degrees(d.omega)
        if 0:
            pylab.figure()
            pylab.plot( dataset.angle + p.zeroAngle , dataset.corr, "-")
            pylab.plot( dom, d.icalc/20000, "o")
            pylab.plot( dom, -d.icalc*d.plr/20000, "o")
            pylab.plot( dom, -d.icalc/d.plr/20000, "o")
            pylab.plot( obs, p.heights,"+")
            for i,j in pairs:
                pylab.plot( [obs[j],dom[i] ],[p.heights[j],d.icalc[i]/20000],"k-")
            pylab.show()
        self.d = d
        self.p = p
        out=open("pairs.dat","w")
        for i,j in pairs:
            out.write("%d %d %d %f %f %f %f %f\n"%(
            d.h[i,0],d.h[i,1],d.h[i,2],dom[i],d.icalc[i],obs[j],p.heights[j],d.plr[j]))
        out.close()
        self.p=p
        self.pairs=pairs

    def gof(self, args, debug = False):
        rx, ry, rz, w, axis = args
        u0 = np.dot( np.dot(  rotmatx(rx),rotmaty(ry)),rotmatz(rz))
        if debug:
            print "orientation",rx,ry,rz,w,axis
            print u0 
        self.d.crystal.set_wvln( w )       
        self.d.axistilt = axis
        self.d.crystal.set_orientation( u0 )
        self.d.crystal.generate_gve(  )
        self.d.computek()
        self.d.computeomegas()
        dom = np.degrees( self.d.omega )
        diff = np.array([dom[i] - self.obspeaks[j] for i,j in self.pairs])
        sum2 = (diff*diff).sum()
        if debug:
#           pylab.plot([dom[i] for i,j in self.pairs], diff,"o")
           pylab.hist(diff,bins=64)
        return np.sqrt(sum2)*100

    def fitrot(self):
        guess = self.start + [ self.d.crystal.wvln, 0.0]
        inc   = [np.radians(0.1),]*3 + [0.001, 0.01]
        pylab.figure()
        print "Before",self.gof(guess, debug=True)
        s = Simplex( self.gof, guess, inc )
        fitted, error, n = s.minimize()    
        print self.gof( fitted ,debug=True )
        pylab.show()
        print np.degrees(fitted[:3]),fitted[3]
    
    
if __name__=="__main__":
    import sys
    energy = float(sys.argv[1])
    ubifile = sys.argv[2]
    peaksfile = sys.argv[3]
    s = simplexfitter(energy, ubifile, peaksfile )
    s.fitrot( )
@}



\section{User interface }

This needs to be fixed. Life is short.
For now we will just hack something together using matplotlib 
and the command line.

@d UIgetrange
@{

def UIgetrange( title, x, ylist):
    lowhigh = []
    pylab.ioff()
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
            ax.set_title("You can close it now")
            fig.canvas.draw()
    cid = fig.canvas.mpl_connect('button_press_event', cb )
    pylab.show()
    return min(lowhigh[-2:]), max(lowhigh[-2:])

@}

Now to pick up a single click.

@d UIgetxy
@{
def UIgetxy( title, x, ylist):
    ret = [0,0]
    pylab.ioff()
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title)
    for y in ylist:
        line, = ax.plot(x, y, ",")
    def cb(e):
        """ call back to record clicks """
        ret[:] = e.xdata, e.ydata
        ax.plot( [e.xdata,e.xdata],ax.get_ylim(),"-")
        ax.plot( ax.get_xlim(),[e.ydata,e.ydata],"-")
        ax.set_title("You can close it now")
        fig.canvas.draw()
    cid = fig.canvas.mpl_connect('button_press_event', cb )
    pylab.show()
    return ret
@}

@d UI 
@{ 
@<UIgetrange@>
@<UIgetxy@>
@}

\subsection{ Running all the programs in order }

@O WWM_process.py
@{

import sys, os, shutil

specfile = sys.argv[1]
dummypars = sys.argv[2]
parfile = specfile + ".json"
pkpars = specfile + "_pks.json"
scans = [int(i) for i in sys.argv[3:]]
print scans

print "Usage: specfile dummypars darknum floodnum scan0 scan1 ..."
print "for example: WWM_process.py wwm_42kev_23092013.spc wwmfrz.json  0 1 2 3 4"


def ex(s):
    print s
    os.system(s)

# copy z:\crystallography\jon\sept13\wwm_edges\wwm_42kev_23092013.spc .
shutil.copyfile("wwmhrz.json",parfile)
ex("python WWMcalib.py %s %s %d %d"%(
    specfile, parfile, scans[0], scans[1]))
s = "python WWMmut.py %s %s %d %d %d"%(
    specfile, parfile, scans[2],scans[-1],(scans[-1]-scans[2])/2+scans[2])
ex(s)
s = "python WWMpeaksearch.py %s %s %s "%(specfile, parfile, pkpars)
s += " ".join([str(i) for i in scans[2:]])
ex(s)
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
        blob_moments, s_1, s_I, mx_I, mn_f, mx_s
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

\section{ Tabulated Fhkls for silicon }

@o si.fhkl
@{
  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

  1  1  1   1.41   2.82   3.1354  0.1017     0.60   1  1  1    1000.0             58.1       -58.10         0.00  180.00  81652.75
  2  2  0   2.30   4.60   1.9201  0.2713     1.61   2  2  0     706.6             65.2       -65.17         0.00  180.00 12 618.85
  3  1  1   2.70   5.40   1.6374  0.3730     2.22   3  1  1     428.5             42.1       -42.10         0.00  180.00 24 449.66
  2  2  2   2.82   5.64   1.5677  0.4069     2.42   2  2  2       0.0              0.0         0.00         0.00    0.00  8 412.07
  4  0  0   3.26   6.51   1.3577  0.5425     3.22   4  0  0     114.3             52.5       -52.48         0.00  180.00  6 308.68
  3  3  1   3.55   7.09   1.2459  0.6442     3.83   3  3  1     167.3             34.6        34.62         0.00    0.00 24 259.70
  4  2  2   3.99   7.98   1.1085  0.8138     4.84   4  2  2     211.6             43.8        43.79         0.00    0.00 24 205.29
  5  1  1   4.23   8.46   1.0451  0.9155     5.44   5  1  1      82.5             29.0        29.01         0.00    0.00 24 182.31
  3  3  3   4.23   8.46   1.0451  0.9155     5.44   3  3  3      27.5             29.0        29.01         0.00    0.00  8 182.31
  4  4  0   4.61   9.21   0.9600  1.0850     6.45   4  4  0      56.1             36.9        36.86         0.00    0.00 12 153.60
  5  3  1   4.82   9.63   0.9180  1.1867     7.05   5  3  1      90.3             24.5        24.47         0.00    0.00 48 140.30
  4  4  2   4.89   9.77   0.9051  1.2206     7.25   4  4  2       0.0              0.0         0.00         0.00    0.00 24 136.36
  6  2  0   5.15  10.30   0.8587  1.3563     8.06   6  2  0      64.1             31.2        31.19         0.00    0.00 24 122.58
  5  3  3   5.34  10.68   0.8282  1.4580     8.66   5  3  3      26.4             20.7       -20.75         0.00  180.00 24 113.93
  6  2  2   5.40  10.81   0.8187  1.4919     8.87   6  2  2       0.0              0.0         0.00         0.00    0.00 24 111.30
  4  4  4   5.64  11.29   0.7839  1.6275     9.67   4  4  4      12.9             26.5       -26.54         0.00  180.00  8 101.90
  5  5  1   5.82  11.64   0.7605  1.7292    10.28   5  5  1      16.1             17.7       -17.68         0.00  180.00 24  95.82
  7  1  1   5.82  11.64   0.7605  1.7292    10.28   7  1  1      16.1             17.7        17.68         0.00    0.00 24  95.82
  6  4  2   6.10  12.20   0.7257  1.8988    11.28   6  4  2      48.2             22.7       -22.69         0.00  180.00 48  87.14
  5  5  3   6.26  12.52   0.7070  2.0005    11.89   5  5  3      10.2             15.1       -15.15         0.00  180.00 24  82.63
  7  3  1   6.26  12.52   0.7070  2.0005    11.89   7  3  1      20.4             15.1       -15.15         0.00  180.00 48  82.63
  8  0  0   6.52  13.04   0.6788  2.1700    12.90   8  0  0       3.9             19.5        19.50         0.00    0.00  6  76.06
  7  3  3   6.67  13.34   0.6635  2.2717    13.50   7  3  3       6.6             13.0       -13.04         0.00  180.00 24  72.59
  6  4  4   6.72  13.44   0.6586  2.3056    13.70   6  4  4       0.0              0.0         0.00         0.00    0.00 24  71.50
  6  6  0   6.92  13.84   0.6400  2.4413    14.51   6  6  0       5.1             16.8       -16.84         0.00  180.00 12  67.45
  8  2  2   6.92  13.84   0.6400  2.4413    14.51   8  2  2      10.3             16.8       -16.84         0.00  180.00 24  67.45
  7  5  1   7.06  14.12   0.6271  2.5430    15.11   7  5  1       8.9             11.3       -11.29         0.00  180.00 48  64.69
  5  5  5   7.06  14.12   0.6271  2.5430    15.11   5  5  5       1.5             11.3        11.29         0.00    0.00  8  64.69
  6  6  2   7.11  14.22   0.6229  2.5769    15.31   6  6  2       0.0              0.0         0.00         0.00    0.00 24  63.82
  8  4  0   7.29  14.59   0.6072  2.7125    16.12   8  4  0       7.0             14.6       -14.62         0.00  180.00 24  60.56
  7  5  3   7.43  14.86   0.5961  2.8142    16.72   7  5  3       6.0              9.8         9.81         0.00    0.00 48  58.32
  9  1  1   7.43  14.86   0.5961  2.8142    16.72   9  1  1       3.0              9.8        -9.81         0.00  180.00 24  58.32
  8  4  2   7.48  14.95   0.5925  2.8481    16.93   8  4  2       0.0              0.0         0.00         0.00    0.00 48  57.60
  6  6  4   7.65  15.30   0.5789  2.9838    17.73   6  6  4       4.8             12.7        12.75         0.00    0.00 24  54.92
  9  3  1   7.78  15.56   0.5693  3.0855    18.34   9  3  1       4.2              8.6        -8.57         0.00  180.00 48  53.06
  8  4  4   7.99  15.99   0.5543  3.2550    19.34   8  4  4       3.4             11.2        11.17         0.00    0.00 24  50.22
  9  3  3   8.12  16.24   0.5458  3.3567    19.95   9  3  3       1.5              7.5         7.53         0.00    0.00 24  48.66
  7  7  1   8.12  16.24   0.5458  3.3567    19.95   7  7  1       1.5              7.5         7.53         0.00    0.00 24  48.66
  7  5  5   8.12  16.24   0.5458  3.3567    19.95   7  5  5       1.5              7.5         7.53         0.00    0.00 24  48.66
 10  2  0   8.32  16.65   0.5325  3.5263    20.96  10  2  0       2.4              9.8        -9.83         0.00  180.00 24  46.25
  8  6  2   8.32  16.65   0.5325  3.5263    20.96   8  6  2       4.8              9.8         9.83         0.00    0.00 48  46.25
  9  5  1   8.44  16.89   0.5250  3.6280    21.56   9  5  1       2.1              6.6         6.63         0.00    0.00 48  44.91
  7  7  3   8.44  16.89   0.5250  3.6280    21.56   7  7  3       1.1              6.6         6.63         0.00    0.00 24  44.91
 10  2  2   8.48  16.97   0.5226  3.6619    21.76  10  2  2       0.0              0.0         0.00         0.00    0.00 24  44.48
  6  6  6   8.48  16.97   0.5226  3.6619    21.76   6  6  6       0.0              0.0         0.00         0.00    0.00  8  44.48
  9  5  3   8.76  17.51   0.5064  3.8992    23.17   9  5  3       1.5              5.9         5.87         0.00    0.00 48  41.69
  8  6  4   8.79  17.59   0.5042  3.9331    23.37   8  6  4       0.0              0.0         0.00         0.00    0.00 48  41.31
 10  4  2   8.95  17.89   0.4958  4.0688    24.18  10  4  2       2.5              7.7         7.71         0.00    0.00 48  39.89
  7  7  5   9.06  18.12   0.4897  4.1705    24.78   7  7  5       0.6              5.2        -5.22         0.00  180.00 24  38.88
 11  1  1   9.06  18.12   0.4897  4.1705    24.78  11  1  1       0.6              5.2        -5.22         0.00  180.00 24  38.88
  8  8  0   9.24  18.48   0.4800  4.3400    25.79   8  8  0       0.5              6.9         6.87         0.00    0.00 12  37.31
  9  5  5   9.35  18.70   0.4745  4.4417    26.40   9  5  5       0.4              4.7        -4.65         0.00  180.00 24  36.42
 11  3  1   9.35  18.70   0.4745  4.4417    26.40  11  3  1       0.8              4.7         4.65         0.00    0.00 48  36.42
  9  7  1   9.35  18.70   0.4745  4.4417    26.40   9  7  1       0.8              4.7         4.65         0.00    0.00 48  36.42
 10  4  4   9.39  18.77   0.4727  4.4756    26.60  10  4  4       0.0              0.0         0.00         0.00    0.00 24  36.13
  8  8  2   9.39  18.77   0.4727  4.4756    26.60   8  8  2       0.0              0.0         0.00         0.00    0.00 24  36.13
 10  6  0   9.53  19.06   0.4657  4.6113    27.40  10  6  0       0.7              6.1         6.14         0.00    0.00 24  35.03
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

  8  6  6   9.53  19.06   0.4657  4.6113    27.40   8  6  6       0.7              6.1        -6.14         0.00  180.00 24  35.03
 11  3  3   9.63  19.27   0.4606  4.7130    28.01  11  3  3       0.3              4.2         4.16         0.00    0.00 24  34.24
  9  7  3   9.63  19.27   0.4606  4.7130    28.01   9  7  3       0.6              4.2        -4.16         0.00  180.00 48  34.24
 10  6  2   9.67  19.34   0.4590  4.7469    28.21  10  6  2       0.0              0.0         0.00         0.00    0.00 48  33.99
 12  0  0   9.81  19.62   0.4526  4.8825    29.02  12  0  0       0.1              5.5        -5.50         0.00  180.00  6  33.00
  8  8  4   9.81  19.62   0.4526  4.8825    29.02   8  8  4       0.5              5.5        -5.50         0.00  180.00 24  33.00
 11  5  1   9.91  19.82   0.4479  4.9842    29.62  11  5  1       0.5              3.7         3.74         0.00    0.00 48  32.30
  7  7  7   9.91  19.82   0.4479  4.9842    29.62   7  7  7       0.1              3.7        -3.74         0.00  180.00  8  32.30
 12  2  2  10.08  20.16   0.4405  5.1538    30.63  12  2  2       0.4              4.9         4.95         0.00    0.00 24  31.19
 10  6  4  10.08  20.16   0.4405  5.1538    30.63  10  6  4       0.8              4.9        -4.95         0.00  180.00 48  31.19
  9  7  5  10.18  20.36   0.4362  5.2555    31.23   9  7  5       0.4              3.4        -3.36         0.00  180.00 48  30.56
 11  5  3  10.18  20.36   0.4362  5.2555    31.23  11  5  3       0.4              3.4        -3.36         0.00  180.00 48  30.56
 12  4  0  10.34  20.69   0.4293  5.4250    32.24  12  4  0       0.3              4.5         4.46         0.00    0.00 24  29.56
  9  9  1  10.44  20.88   0.4254  5.5267    32.84   9  9  1       0.1              3.0        -3.03         0.00  180.00 24  28.99
  8  8  6  10.47  20.95   0.4241  5.5606    33.05   8  8  6       0.0              0.0         0.00         0.00    0.00 24  28.81
 12  4  2  10.47  20.95   0.4241  5.5606    33.05  12  4  2       0.0              0.0         0.00         0.00    0.00 48  28.81
 10  8  2  10.60  21.20   0.4190  5.6963    33.85  10  8  2       0.5              4.0        -4.03         0.00  180.00 48  28.09
 11  7  1  10.70  21.40   0.4153  5.7980    34.46  11  7  1       0.2              2.7        -2.75         0.00  180.00 48  27.57
 11  5  5  10.70  21.40   0.4153  5.7980    34.46  11  5  5       0.1              2.7        -2.75         0.00  180.00 24  27.57
  9  9  3  10.70  21.40   0.4153  5.7980    34.46   9  9  3       0.1              2.7        -2.75         0.00  180.00 24  27.57
 13  1  1  10.70  21.40   0.4153  5.7980    34.46  13  1  1       0.1              2.7         2.75         0.00    0.00 24  27.57
 10  6  6  10.73  21.46   0.4141  5.8319    34.66  10  6  6       0.0              0.0         0.00         0.00    0.00 24  27.40
 12  4  4  10.85  21.71   0.4094  5.9675    35.46  12  4  4       0.2              3.7        -3.65         0.00  180.00 24  26.75
 11  7  3  10.95  21.90   0.4059  6.0692    36.07  11  7  3       0.2              2.5        -2.49         0.00  180.00 48  26.28
 13  3  1  10.95  21.90   0.4059  6.0692    36.07  13  3  1       0.2              2.5         2.49         0.00    0.00 48  26.28
  9  7  7  10.95  21.90   0.4059  6.0692    36.07   9  7  7       0.1              2.5         2.49         0.00    0.00 24  26.28
 10  8  4  10.98  21.96   0.4048  6.1031    36.27  10  8  4       0.0              0.0         0.00         0.00    0.00 48  26.12
 12  6  2  11.10  22.20   0.4004  6.2388    37.08  12  6  2       0.3              3.3        -3.31         0.00  180.00 48  25.52
 13  3  3  11.19  22.39   0.3971  6.3405    37.68  13  3  3       0.1              2.3        -2.26         0.00  180.00 24  25.09
  9  9  5  11.19  22.39   0.3971  6.3405    37.68   9  9  5       0.1              2.3         2.26         0.00    0.00 24  25.09
  8  8  8  11.34  22.69   0.3919  6.5100    38.69   8  8  8       0.0              3.0         3.02         0.00    0.00  8  24.40
 13  5  1  11.43  22.87   0.3889  6.6117    39.29  13  5  1       0.1              2.1        -2.06         0.00  180.00 48  24.00
 11  7  5  11.43  22.87   0.3889  6.6117    39.29  11  7  5       0.1              2.1         2.06         0.00    0.00 48  24.00
 12  6  4  11.46  22.93   0.3879  6.6456    39.49  12  6  4       0.0              0.0         0.00         0.00    0.00 48  23.88
 14  2  0  11.58  23.16   0.3840  6.7813    40.30  14  2  0       0.1              2.7         2.75         0.00    0.00 24  23.37
 10 10  0  11.58  23.16   0.3840  6.7813    40.30  10 10  0       0.0              2.7        -2.75         0.00  180.00 12  23.37
 10  8  6  11.58  23.16   0.3840  6.7813    40.30  10  8  6       0.2              2.7         2.75         0.00    0.00 48  23.37
 13  5  3  11.67  23.34   0.3812  6.8830    40.90  13  5  3       0.1              1.9        -1.88         0.00  180.00 48  23.00
 11  9  1  11.67  23.34   0.3812  6.8830    40.90  11  9  1       0.1              1.9        -1.88         0.00  180.00 48  23.00
 14  2  2  11.70  23.40   0.3802  6.9169    41.11  14  2  2       0.0              0.0         0.00         0.00    0.00 24  22.88
 10 10  2  11.70  23.40   0.3802  6.9169    41.11  10 10  2       0.0              0.0         0.00         0.00    0.00 24  22.88
 12  8  0  11.81  23.63   0.3766  7.0525    41.91  12  8  0       0.1              2.5        -2.51         0.00  180.00 24  22.42
  9  9  7  11.90  23.80   0.3739  7.1542    42.52   9  9  7       0.0              1.7         1.71         0.00    0.00 24  22.08
 11  9  3  11.90  23.80   0.3739  7.1542    42.52  11  9  3       0.1              1.7         1.71         0.00    0.00 48  22.08
 12  8  2  11.93  23.86   0.3730  7.1881    42.72  12  8  2       0.0              0.0         0.00         0.00    0.00 48  21.97
 10 10  4  12.04  24.08   0.3695  7.3238    43.52  10 10  4       0.1              2.3         2.29         0.00    0.00 24  21.54
 14  4  2  12.04  24.08   0.3695  7.3238    43.52  14  4  2       0.1              2.3        -2.29         0.00  180.00 48  21.54
 12  6  6  12.04  24.08   0.3695  7.3238    43.52  12  6  6       0.1              2.3         2.29         0.00    0.00 24  21.54
 11  7  7  12.13  24.25   0.3670  7.4255    44.13  11  7  7       0.0              1.6         1.57         0.00    0.00 24  21.22
 13  7  1  12.13  24.25   0.3670  7.4255    44.13  13  7  1       0.1              1.6        -1.57         0.00  180.00 48  21.22
 13  5  5  12.13  24.25   0.3670  7.4255    44.13  13  5  5       0.0              1.6         1.57         0.00    0.00 24  21.22
 12  8  4  12.27  24.53   0.3629  7.5950    45.14  12  8  4       0.1              2.1         2.10         0.00    0.00 48  20.72
 11  9  5  12.35  24.70   0.3605  7.6967    45.74  11  9  5       0.0              1.4         1.43         0.00    0.00 48  20.43
 15  1  1  12.35  24.70   0.3605  7.6967    45.74  15  1  1       0.0              1.4         1.43         0.00    0.00 24  20.43
 13  7  3  12.35  24.70   0.3605  7.6967    45.74  13  7  3       0.0              1.4         1.43         0.00    0.00 48  20.43
 14  4  4  12.38  24.75   0.3597  7.7306    45.94  14  4  4       0.0              0.0         0.00         0.00    0.00 24  20.33
 10  8  8  12.38  24.75   0.3597  7.7306    45.94  10  8  8       0.0              0.0         0.00         0.00    0.00 24  20.33
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 14  6  0  12.49  24.97   0.3565  7.8663    46.75  14  6  0       0.0              1.9        -1.92         0.00  180.00 24  19.96
 15  3  1  12.57  25.14   0.3543  7.9680    47.35  15  3  1       0.0              1.3        -1.32         0.00  180.00 48  19.68
 14  6  2  12.60  25.19   0.3535  8.0019    47.55  14  6  2       0.0              0.0         0.00         0.00    0.00 48  19.60
 10 10  6  12.60  25.19   0.3535  8.0019    47.55  10 10  6       0.0              0.0         0.00         0.00    0.00 24  19.60
 15  3  3  12.78  25.57   0.3484  8.2392    48.96  15  3  3       0.0              1.2        -1.21         0.00  180.00 24  18.99
 11 11  1  12.78  25.57   0.3484  8.2392    48.96  11 11  1       0.0              1.2         1.21         0.00    0.00 24  18.99
 13  7  5  12.78  25.57   0.3484  8.2392    48.96  13  7  5       0.0              1.2         1.21         0.00    0.00 48  18.99
  9  9  9  12.78  25.57   0.3484  8.2392    48.96   9  9  9       0.0              1.2        -1.21         0.00  180.00  8  18.99
 12  8  6  12.81  25.62   0.3477  8.2731    49.17  12  8  6       0.0              0.0         0.00         0.00    0.00 48  18.91
 12 10  2  12.92  25.83   0.3449  8.4088    49.97  12 10  2       0.1              1.6         1.62         0.00    0.00 48  18.58
 14  6  4  12.92  25.83   0.3449  8.4088    49.97  14  6  4       0.1              1.6         1.62         0.00    0.00 48  18.58
 11  9  7  13.00  25.99   0.3428  8.5105    50.58  11  9  7       0.0              1.1        -1.11         0.00  180.00 48  18.34
 11 11  3  13.00  25.99   0.3428  8.5105    50.58  11 11  3       0.0              1.1         1.11         0.00    0.00 24  18.34
 15  5  1  13.00  25.99   0.3428  8.5105    50.58  15  5  1       0.0              1.1        -1.11         0.00  180.00 48  18.34
 13  9  1  13.00  25.99   0.3428  8.5105    50.58  13  9  1       0.0              1.1         1.11         0.00    0.00 48  18.34
 16  0  0  13.13  26.26   0.3394  8.6800    51.58  16  0  0       0.0              1.5         1.49         0.00    0.00  6  17.96
 13  9  3  13.21  26.41   0.3375  8.7817    52.19  13  9  3       0.0              1.0         1.02         0.00    0.00 48  17.73
 15  5  3  13.21  26.41   0.3375  8.7817    52.19  15  5  3       0.0              1.0         1.02         0.00    0.00 48  17.73
 12 10  4  13.23  26.46   0.3368  8.8156    52.39  12 10  4       0.0              0.0         0.00         0.00    0.00 48  17.66
 16  2  2  13.33  26.67   0.3342  8.9513    53.20  16  2  2       0.0              1.4        -1.37         0.00  180.00 24  17.37
 10 10  8  13.33  26.67   0.3342  8.9513    53.20  10 10  8       0.0              1.4        -1.37         0.00  180.00 24  17.37
 14  8  2  13.33  26.67   0.3342  8.9513    53.20  14  8  2       0.0              1.4         1.37         0.00    0.00 48  17.37
 11 11  5  13.41  26.82   0.3324  9.0530    53.80  11 11  5       0.0              0.9        -0.94         0.00  180.00 24  17.16
 13  7  7  13.41  26.82   0.3324  9.0530    53.80  13  7  7       0.0              0.9        -0.94         0.00  180.00 24  17.16
 14  6  6  13.44  26.87   0.3317  9.0869    54.00  14  6  6       0.0              0.0         0.00         0.00    0.00 24  17.09
 16  4  0  13.54  27.08   0.3293  9.2225    54.81  16  4  0       0.0              1.3        -1.26         0.00  180.00 24  16.82
 12  8  8  13.54  27.08   0.3293  9.2225    54.81  12  8  8       0.0              1.3        -1.26         0.00  180.00 24  16.82
 15  7  1  13.62  27.23   0.3275  9.3242    55.41  15  7  1       0.0              0.9         0.86         0.00    0.00 48  16.62
 15  5  5  13.62  27.23   0.3275  9.3242    55.41  15  5  5       0.0              0.9         0.86         0.00    0.00 24  16.62
 13  9  5  13.62  27.23   0.3275  9.3242    55.41  13  9  5       0.0              0.9        -0.86         0.00  180.00 48  16.62
 16  4  2  13.64  27.28   0.3269  9.3581    55.61  16  4  2       0.0              0.0         0.00         0.00    0.00 48  16.56
 14  8  4  13.64  27.28   0.3269  9.3581    55.61  14  8  4       0.0              0.0         0.00         0.00    0.00 48  16.56
 12 10  6  13.74  27.48   0.3245  9.4938    56.42  12 10  6       0.0              1.2        -1.16         0.00  180.00 48  16.30
 15  7  3  13.82  27.63   0.3228  9.5955    57.02  15  7  3       0.0              0.8         0.80         0.00    0.00 48  16.12
 11  9  9  13.82  27.63   0.3228  9.5955    57.02  11  9  9       0.0              0.8        -0.80         0.00  180.00 24  16.12
 16  4  4  13.94  27.88   0.3200  9.7650    58.03  16  4  4       0.0              1.1         1.07         0.00    0.00 24  15.81
 12 12  0  13.94  27.88   0.3200  9.7650    58.03  12 12  0       0.0              1.1         1.07         0.00    0.00 12  15.81
 17  1  1  14.01  28.03   0.3184  9.8667    58.64  17  1  1       0.0              0.7        -0.73         0.00  180.00 24  15.64
 13 11  1  14.01  28.03   0.3184  9.8667    58.64  13 11  1       0.0              0.7         0.73         0.00    0.00 48  15.64
 11 11  7  14.01  28.03   0.3184  9.8667    58.64  11 11  7       0.0              0.7        -0.73         0.00  180.00 24  15.64
 12 12  2  14.04  28.08   0.3178  9.9006    58.84  12 12  2       0.0              0.0         0.00         0.00    0.00 24  15.58
 16  6  2  14.14  28.27   0.3157 10.0363    59.64  16  6  2       0.0              1.0         0.99         0.00    0.00 48  15.35
 14 10  0  14.14  28.27   0.3157 10.0363    59.64  14 10  0       0.0              1.0         0.99         0.00    0.00 24  15.35
 14  8  6  14.14  28.27   0.3157 10.0363    59.64  14  8  6       0.0              1.0        -0.99         0.00  180.00 48  15.35
 13  9  7  14.21  28.42   0.3141 10.1380    60.25  13  9  7       0.0              0.7        -0.68         0.00  180.00 48  15.18
 17  3  1  14.21  28.42   0.3141 10.1380    60.25  17  3  1       0.0              0.7        -0.68         0.00  180.00 48  15.18
 13 11  3  14.21  28.42   0.3141 10.1380    60.25  13 11  3       0.0              0.7        -0.68         0.00  180.00 48  15.18
 15  7  5  14.21  28.42   0.3141 10.1380    60.25  15  7  5       0.0              0.7        -0.68         0.00  180.00 48  15.18
 14 10  2  14.23  28.47   0.3135 10.1719    60.45  14 10  2       0.0              0.0         0.00         0.00    0.00 48  15.13
 10 10 10  14.23  28.47   0.3135 10.1719    60.45  10 10 10       0.0              0.0         0.00         0.00    0.00  8  15.13
 12 12  4  14.33  28.66   0.3115 10.3075    61.26  12 12  4       0.0              0.9        -0.91         0.00  180.00 24  14.91
 15  9  1  14.40  28.80   0.3099 10.4092    61.86  15  9  1       0.0              0.6         0.63         0.00    0.00 48  14.75
 17  3  3  14.40  28.80   0.3099 10.4092    61.86  17  3  3       0.0              0.6         0.63         0.00    0.00 24  14.75
 12 10  8  14.43  28.85   0.3094 10.4431    62.06  12 10  8       0.0              0.0         0.00         0.00    0.00 48  14.70
 16  6  4  14.43  28.85   0.3094 10.4431    62.06  16  6  4       0.0              0.0         0.00         0.00    0.00 48  14.70
 14 10  4  14.52  29.04   0.3075 10.5788    62.87  14 10  4       0.0              0.8        -0.84         0.00  180.00 48  14.50
 17  5  1  14.59  29.18   0.3060 10.6805    63.47  17  5  1       0.0              0.6         0.58         0.00    0.00 48  14.34
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 15  9  3  14.59  29.18   0.3060 10.6805    63.47  15  9  3       0.0              0.6        -0.58         0.00  180.00 48  14.34
 13 11  5  14.59  29.18   0.3060 10.6805    63.47  13 11  5       0.0              0.6        -0.58         0.00  180.00 48  14.34
 16  8  0  14.71  29.42   0.3036 10.8500    64.48  16  8  0       0.0              0.8         0.78         0.00    0.00 24  14.10
 17  5  3  14.78  29.56   0.3022 10.9517    65.08  17  5  3       0.0              0.5         0.53         0.00    0.00 48  13.96
 11 11  9  14.78  29.56   0.3022 10.9517    65.08  11 11  9       0.0              0.5         0.53         0.00    0.00 24  13.96
 15  7  7  14.78  29.56   0.3022 10.9517    65.08  15  7  7       0.0              0.5        -0.53         0.00  180.00 24  13.96
 16  8  2  14.80  29.61   0.3017 10.9856    65.29  16  8  2       0.0              0.0         0.00         0.00    0.00 48  13.91
 12 12  6  14.80  29.61   0.3017 10.9856    65.29  12 12  6       0.0              0.0         0.00         0.00    0.00 24  13.91
 14  8  8  14.80  29.61   0.3017 10.9856    65.29  14  8  8       0.0              0.0         0.00         0.00    0.00 24  13.91
 18  2  0  14.90  29.79   0.2999 11.1213    66.09  18  2  0       0.0              0.7        -0.72         0.00  180.00 24  13.72
 16  6  6  14.90  29.79   0.2999 11.1213    66.09  16  6  6       0.0              0.7        -0.72         0.00  180.00 24  13.72
 15  9  5  14.97  29.93   0.2985 11.2230    66.70  15  9  5       0.0              0.5        -0.49         0.00  180.00 48  13.59
 13  9  9  14.97  29.93   0.2985 11.2230    66.70  13  9  9       0.0              0.5         0.49         0.00    0.00 24  13.59
 14 10  6  14.99  29.98   0.2981 11.2569    66.90  14 10  6       0.0              0.0         0.00         0.00    0.00 48  13.54
 18  2  2  14.99  29.98   0.2981 11.2569    66.90  18  2  2       0.0              0.0         0.00         0.00    0.00 24  13.54
 16  8  4  15.08  30.16   0.2963 11.3925    67.70  16  8  4       0.0              0.7        -0.67         0.00  180.00 48  13.37
 13 11  7  15.15  30.30   0.2950 11.4942    68.31  13 11  7       0.0              0.5         0.46         0.00    0.00 48  13.24
 17  7  1  15.15  30.30   0.2950 11.4942    68.31  17  7  1       0.0              0.5         0.46         0.00    0.00 48  13.24
 17  5  5  15.15  30.30   0.2950 11.4942    68.31  17  5  5       0.0              0.5        -0.46         0.00  180.00 24  13.24
 13 13  1  15.15  30.30   0.2950 11.4942    68.31  13 13  1       0.0              0.5        -0.46         0.00  180.00 24  13.24
 14 12  2  15.26  30.53   0.2928 11.6638    69.32  14 12  2       0.0              0.6        -0.62         0.00  180.00 48  13.02
 18  4  2  15.26  30.53   0.2928 11.6638    69.32  18  4  2       0.0              0.6         0.62         0.00    0.00 48  13.02
 12 10 10  15.26  30.53   0.2928 11.6638    69.32  12 10 10       0.0              0.6         0.62         0.00    0.00 24  13.02
 17  7  3  15.33  30.67   0.2915 11.7655    69.92  17  7  3       0.0              0.4        -0.42         0.00  180.00 48  12.90
 13 13  3  15.33  30.67   0.2915 11.7655    69.92  13 13  3       0.0              0.4        -0.42         0.00  180.00 24  12.90
 15 11  1  15.33  30.67   0.2915 11.7655    69.92  15 11  1       0.0              0.4        -0.42         0.00  180.00 48  12.90
 12 12  8  15.45  30.89   0.2895 11.9350    70.93  12 12  8       0.0              0.6         0.57         0.00    0.00 24  12.70
 15  9  7  15.51  31.03   0.2882 12.0367    71.53  15  9  7       0.0              0.4         0.39         0.00    0.00 48  12.58
 15 11  3  15.51  31.03   0.2882 12.0367    71.53  15 11  3       0.0              0.4        -0.39         0.00  180.00 48  12.58
 14 12  4  15.54  31.07   0.2878 12.0706    71.73  14 12  4       0.0              0.0         0.00         0.00    0.00 48  12.54
 16  8  6  15.54  31.07   0.2878 12.0706    71.73  16  8  6       0.0              0.0         0.00         0.00    0.00 48  12.54
 18  4  4  15.54  31.07   0.2878 12.0706    71.73  18  4  4       0.0              0.0         0.00         0.00    0.00 24  12.54
 18  6  0  15.62  31.25   0.2862 12.2063    72.54  18  6  0       0.0              0.5         0.53         0.00    0.00 24  12.39
 14 10  8  15.62  31.25   0.2862 12.2063    72.54  14 10  8       0.0              0.5         0.53         0.00    0.00 48  12.39
 16 10  2  15.62  31.25   0.2862 12.2063    72.54  16 10  2       0.0              0.5        -0.53         0.00  180.00 48  12.39
 17  7  5  15.69  31.38   0.2850 12.3080    73.14  17  7  5       0.0              0.4        -0.36         0.00  180.00 48  12.28
 13 13  5  15.69  31.38   0.2850 12.3080    73.14  13 13  5       0.0              0.4         0.36         0.00    0.00 24  12.28
 19  1  1  15.69  31.38   0.2850 12.3080    73.14  19  1  1       0.0              0.4        -0.36         0.00  180.00 24  12.28
 11 11 11  15.69  31.38   0.2850 12.3080    73.14  11 11 11       0.0              0.4         0.36         0.00    0.00  8  12.28
 18  6  2  15.71  31.43   0.2846 12.3419    73.35  18  6  2       0.0              0.0         0.00         0.00    0.00 48  12.24
 13 11  9  15.87  31.74   0.2820 12.5792    74.76  13 11  9       0.0              0.3         0.34         0.00    0.00 48  11.98
 17  9  1  15.87  31.74   0.2820 12.5792    74.76  17  9  1       0.0              0.3        -0.34         0.00  180.00 48  11.98
 19  3  1  15.87  31.74   0.2820 12.5792    74.76  19  3  1       0.0              0.3         0.34         0.00    0.00 48  11.98
 15 11  5  15.87  31.74   0.2820 12.5792    74.76  15 11  5       0.0              0.3         0.34         0.00    0.00 48  11.98
 16 10  4  15.89  31.78   0.2816 12.6131    74.96  16 10  4       0.0              0.0         0.00         0.00    0.00 48  11.95
 14 12  6  15.98  31.95   0.2801 12.7488    75.76  14 12  6       0.0              0.5         0.45         0.00    0.00 48  11.81
 18  6  4  15.98  31.95   0.2801 12.7488    75.76  18  6  4       0.0              0.5        -0.45         0.00  180.00 48  11.81
 19  3  3  16.04  32.08   0.2790 12.8505    76.37  19  3  3       0.0              0.3         0.31         0.00    0.00 24  11.70
 17  9  3  16.04  32.08   0.2790 12.8505    76.37  17  9  3       0.0              0.3        -0.31         0.00  180.00 48  11.70
 16  8  8  16.15  32.30   0.2771 13.0200    77.38  16  8  8       0.0              0.4         0.42         0.00    0.00 24  11.53
 17  7  7  16.22  32.43   0.2761 13.1217    77.98  17  7  7       0.0              0.3         0.29         0.00    0.00 24  11.43
 13 13  7  16.22  32.43   0.2761 13.1217    77.98  13 13  7       0.0              0.3         0.29         0.00    0.00 24  11.43
 15  9  9  16.22  32.43   0.2761 13.1217    77.98  15  9  9       0.0              0.3         0.29         0.00    0.00 24  11.43
 19  5  1  16.22  32.43   0.2761 13.1217    77.98  19  5  1       0.0              0.3         0.29         0.00    0.00 48  11.43
 12 12 10  16.24  32.47   0.2757 13.1556    78.18  12 12 10       0.0              0.0         0.00         0.00    0.00 24  11.40
 18  8  2  16.32  32.65   0.2743 13.2913    78.99  18  8  2       0.0              0.4        -0.39         0.00  180.00 48  11.27
 16 10  6  16.32  32.65   0.2743 13.2913    78.99  16 10  6       0.0              0.4         0.39         0.00    0.00 48  11.27
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 14 14  0  16.32  32.65   0.2743 13.2913    78.99  14 14  0       0.0              0.4        -0.39         0.00  180.00 12  11.27
 15 13  1  16.39  32.77   0.2733 13.3930    79.59  15 13  1       0.0              0.3        -0.27         0.00  180.00 48  11.18
 19  5  3  16.39  32.77   0.2733 13.3930    79.59  19  5  3       0.0              0.3        -0.27         0.00  180.00 48  11.18
 17  9  5  16.39  32.77   0.2733 13.3930    79.59  17  9  5       0.0              0.3         0.27         0.00    0.00 48  11.18
 15 11  7  16.39  32.77   0.2733 13.3930    79.59  15 11  7       0.0              0.3         0.27         0.00    0.00 48  11.18
 14 10 10  16.41  32.82   0.2729 13.4269    79.79  14 10 10       0.0              0.0         0.00         0.00    0.00 24  11.15
 14 14  2  16.41  32.82   0.2729 13.4269    79.79  14 14  2       0.0              0.0         0.00         0.00    0.00 24  11.15
 18  6  6  16.41  32.82   0.2729 13.4269    79.79  18  6  6       0.0              0.0         0.00         0.00    0.00 24  11.15
 16 12  0  16.49  32.99   0.2715 13.5625    80.60  16 12  0       0.0              0.4        -0.36         0.00  180.00 24  11.02
 20  0  0  16.49  32.99   0.2715 13.5625    80.60  20  0  0       0.0              0.4        -0.36         0.00  180.00  6  11.02
 15 13  3  16.56  33.11   0.2705 13.6642    81.20  15 13  3       0.0              0.2         0.25         0.00    0.00 48  10.93
 14 12  8  16.58  33.16   0.2702 13.6981    81.41  14 12  8       0.0              0.0         0.00         0.00    0.00 48  10.90
 18  8  4  16.58  33.16   0.2702 13.6981    81.41  18  8  4       0.0              0.0         0.00         0.00    0.00 48  10.90
 16 12  2  16.58  33.16   0.2702 13.6981    81.41  16 12  2       0.0              0.0         0.00         0.00    0.00 48  10.90
 14 14  4  16.66  33.32   0.2689 13.8338    82.21  14 14  4       0.0              0.3         0.34         0.00    0.00 24  10.78
 20  2  2  16.66  33.32   0.2689 13.8338    82.21  20  2  2       0.0              0.3         0.34         0.00    0.00 24  10.78
 19  5  5  16.73  33.45   0.2679 13.9355    82.82  19  5  5       0.0              0.2        -0.23         0.00  180.00 24  10.69
 17 11  1  16.73  33.45   0.2679 13.9355    82.82  17 11  1       0.0              0.2        -0.23         0.00  180.00 48  10.69
 19  7  1  16.73  33.45   0.2679 13.9355    82.82  19  7  1       0.0              0.2        -0.23         0.00  180.00 48  10.69
 13 11 11  16.73  33.45   0.2679 13.9355    82.82  13 11 11       0.0              0.2        -0.23         0.00  180.00 24  10.69
 20  4  0  16.83  33.66   0.2663 14.1050    83.82  20  4  0       0.0              0.3         0.31         0.00    0.00 24  10.55
 16 12  4  16.83  33.66   0.2663 14.1050    83.82  16 12  4       0.0              0.3         0.31         0.00    0.00 48  10.55
 15 13  5  16.89  33.78   0.2653 14.2067    84.43  15 13  5       0.0              0.2         0.21         0.00    0.00 48  10.46
 19  7  3  16.89  33.78   0.2653 14.2067    84.43  19  7  3       0.0              0.2        -0.21         0.00  180.00 48  10.46
 17  9  7  16.89  33.78   0.2653 14.2067    84.43  17  9  7       0.0              0.2         0.21         0.00    0.00 48  10.46
 13 13  9  16.89  33.78   0.2653 14.2067    84.43  13 13  9       0.0              0.2        -0.21         0.00  180.00 24  10.46
 17 11  3  16.89  33.78   0.2653 14.2067    84.43  17 11  3       0.0              0.2         0.21         0.00    0.00 48  10.46
 16 10  8  16.91  33.83   0.2650 14.2406    84.63  16 10  8       0.0              0.0         0.00         0.00    0.00 48  10.44
 20  4  2  16.91  33.83   0.2650 14.2406    84.63  20  4  2       0.0              0.0         0.00         0.00    0.00 48  10.44
 18  8  6  17.00  33.99   0.2637 14.3763    85.44  18  8  6       0.0              0.3         0.29         0.00    0.00 48  10.33
 18 10  0  17.00  33.99   0.2637 14.3763    85.44  18 10  0       0.0              0.3        -0.29         0.00  180.00 24  10.33
 15 11  9  17.06  34.11   0.2628 14.4780    86.04  15 11  9       0.0              0.2        -0.20         0.00  180.00 48  10.25
 18 10  2  17.08  34.16   0.2625 14.5119    86.24  18 10  2       0.0              0.0         0.00         0.00    0.00 48  10.22
 14 14  6  17.08  34.16   0.2625 14.5119    86.24  14 14  6       0.0              0.0         0.00         0.00    0.00 24  10.22
 20  4  4  17.16  34.32   0.2613 14.6475    87.05  20  4  4       0.0              0.3        -0.27         0.00  180.00 24  10.11
 12 12 12  17.16  34.32   0.2613 14.6475    87.05  12 12 12       0.0              0.3        -0.27         0.00  180.00  8  10.11
 17 11  5  17.22  34.44   0.2604 14.7492    87.65  17 11  5       0.0              0.2         0.18         0.00    0.00 48  10.03
 19  7  5  17.22  34.44   0.2604 14.7492    87.65  19  7  5       0.0              0.2         0.18         0.00    0.00 48  10.03
 16 12  6  17.24  34.48   0.2601 14.7831    87.85  16 12  6       0.0              0.0         0.00         0.00    0.00 48  10.01
 18 10  4  17.32  34.65   0.2589 14.9188    88.66  18 10  4       0.0              0.2         0.25         0.00    0.00 48   9.91
 20  6  2  17.32  34.65   0.2589 14.9188    88.66  20  6  2       0.0              0.2        -0.25         0.00  180.00 48   9.91
 14 12 10  17.32  34.65   0.2589 14.9188    88.66  14 12 10       0.0              0.2        -0.25         0.00  180.00 48   9.91
 21  1  1  17.38  34.77   0.2580 15.0205    89.26  21  1  1       0.0              0.2         0.17         0.00    0.00 24   9.83
 19  9  1  17.38  34.77   0.2580 15.0205    89.26  19  9  1       0.0              0.2        -0.17         0.00  180.00 48   9.83
 15 13  7  17.38  34.77   0.2580 15.0205    89.26  15 13  7       0.0              0.2        -0.17         0.00  180.00 48   9.83
 21  3  1  17.55  35.09   0.2557 15.2917    90.88  21  3  1       0.0              0.2         0.16         0.00    0.00 48   9.63
 19  9  3  17.55  35.09   0.2557 15.2917    90.88  19  9  3       0.0              0.2         0.16         0.00    0.00 48   9.63
 15 15  1  17.55  35.09   0.2557 15.2917    90.88  15 15  1       0.0              0.2         0.16         0.00    0.00 24   9.63
 17  9  9  17.55  35.09   0.2557 15.2917    90.88  17  9  9       0.0              0.2        -0.16         0.00  180.00 24   9.63
 18  8  8  17.57  35.13   0.2554 15.3256    91.08  18  8  8       0.0              0.0         0.00         0.00    0.00 24   9.61
 20  6  4  17.57  35.13   0.2554 15.3256    91.08  20  6  4       0.0              0.0         0.00         0.00    0.00 48   9.61
 16 10 10  17.65  35.29   0.2543 15.4613    91.88  16 10 10       0.0              0.2        -0.22         0.00  180.00 24   9.51
 16 14  2  17.65  35.29   0.2543 15.4613    91.88  16 14  2       0.0              0.2         0.22         0.00    0.00 48   9.51
 14 14  8  17.65  35.29   0.2543 15.4613    91.88  14 14  8       0.0              0.2        -0.22         0.00  180.00 24   9.51
 19  7  7  17.71  35.41   0.2535 15.5630    92.49  19  7  7       0.0              0.1         0.15         0.00    0.00 24   9.44
 17 11  7  17.71  35.41   0.2535 15.5630    92.49  17 11  7       0.0              0.1        -0.15         0.00  180.00 48   9.44
 17 13  1  17.71  35.41   0.2535 15.5630    92.49  17 13  1       0.0              0.1         0.15         0.00    0.00 48   9.44
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 15 15  3  17.71  35.41   0.2535 15.5630    92.49  15 15  3       0.0              0.1         0.15         0.00    0.00 24   9.44
 13 13 11  17.71  35.41   0.2535 15.5630    92.49  13 13 11       0.0              0.1        -0.15         0.00  180.00 24   9.44
 21  3  3  17.71  35.41   0.2535 15.5630    92.49  21  3  3       0.0              0.1        -0.15         0.00  180.00 24   9.44
 18 10  6  17.73  35.45   0.2532 15.5969    92.69  18 10  6       0.0              0.0         0.00         0.00    0.00 48   9.42
 20  8  0  17.80  35.61   0.2521 15.7325    93.50  20  8  0       0.0              0.2        -0.20         0.00  180.00 24   9.33
 16 12  8  17.80  35.61   0.2521 15.7325    93.50  16 12  8       0.0              0.2        -0.20         0.00  180.00 48   9.33
 21  5  1  17.86  35.73   0.2513 15.8342    94.10  21  5  1       0.0              0.1        -0.14         0.00  180.00 48   9.26
 15 11 11  17.86  35.73   0.2513 15.8342    94.10  15 11 11       0.0              0.1        -0.14         0.00  180.00 24   9.26
 19  9  5  17.86  35.73   0.2513 15.8342    94.10  19  9  5       0.0              0.1         0.14         0.00    0.00 48   9.26
 17 13  3  17.86  35.73   0.2513 15.8342    94.10  17 13  3       0.0              0.1         0.14         0.00    0.00 48   9.26
 16 14  4  17.88  35.77   0.2510 15.8681    94.30  16 14  4       0.0              0.0         0.00         0.00    0.00 48   9.24
 20  8  2  17.88  35.77   0.2510 15.8681    94.30  20  8  2       0.0              0.0         0.00         0.00    0.00 48   9.24
 20  6  6  17.96  35.93   0.2500 16.0038    95.11  20  6  6       0.0              0.2         0.19         0.00    0.00 24   9.15
 18 12  2  17.96  35.93   0.2500 16.0038    95.11  18 12  2       0.0              0.2         0.19         0.00    0.00 48   9.15
 15 13  9  18.02  36.04   0.2492 16.1055    95.71  15 13  9       0.0              0.1        -0.13         0.00  180.00 48   9.09
 21  5  3  18.02  36.04   0.2492 16.1055    95.71  21  5  3       0.0              0.1        -0.13         0.00  180.00 48   9.09
 15 15  5  18.02  36.04   0.2492 16.1055    95.71  15 15  5       0.0              0.1        -0.13         0.00  180.00 24   9.09
 20  8  4  18.12  36.24   0.2479 16.2750    96.72  20  8  4       0.0              0.2         0.17         0.00    0.00 48   8.98
 17 13  5  18.18  36.36   0.2471 16.3767    97.32  17 13  5       0.0              0.1        -0.12         0.00  180.00 48   8.91
 19 11  1  18.18  36.36   0.2471 16.3767    97.32  19 11  1       0.0              0.1         0.12         0.00    0.00 48   8.91
 18 12  4  18.20  36.39   0.2469 16.4106    97.53  18 12  4       0.0              0.0         0.00         0.00    0.00 48   8.89
 14 12 12  18.20  36.39   0.2469 16.4106    97.53  14 12 12       0.0              0.0         0.00         0.00    0.00 24   8.89
 18 10  8  18.28  36.55   0.2458 16.5463    98.33  18 10  8       0.0              0.2        -0.16         0.00  180.00 48   8.81
 22  2  0  18.28  36.55   0.2458 16.5463    98.33  22  2  0       0.0              0.2         0.16         0.00    0.00 24   8.81
 16 14  6  18.28  36.55   0.2458 16.5463    98.33  16 14  6       0.0              0.2        -0.16         0.00  180.00 48   8.81
 17 11  9  18.33  36.67   0.2451 16.6480    98.94  17 11  9       0.0              0.1        -0.11         0.00  180.00 48   8.75
 21  7  1  18.33  36.67   0.2451 16.6480    98.94  21  7  1       0.0              0.1        -0.11         0.00  180.00 48   8.75
 21  5  5  18.33  36.67   0.2451 16.6480    98.94  21  5  5       0.0              0.1         0.11         0.00    0.00 24   8.75
 19  9  7  18.33  36.67   0.2451 16.6480    98.94  19  9  7       0.0              0.1        -0.11         0.00  180.00 48   8.75
 19 11  3  18.33  36.67   0.2451 16.6480    98.94  19 11  3       0.0              0.1         0.11         0.00    0.00 48   8.75
 22  2  2  18.35  36.71   0.2448 16.6819    99.14  22  2  2       0.0              0.0         0.00         0.00    0.00 24   8.73
 14 14 10  18.35  36.71   0.2448 16.6819    99.14  14 14 10       0.0              0.0         0.00         0.00    0.00 24   8.73
 15 15  7  18.49  36.97   0.2431 16.9192   100.55  15 15  7       0.0              0.1        -0.10         0.00  180.00 24   8.59
 21  7  3  18.49  36.97   0.2431 16.9192   100.55  21  7  3       0.0              0.1         0.10         0.00    0.00 48   8.59
 16 12 10  18.51  37.01   0.2429 16.9531   100.75  16 12 10       0.0              0.0         0.00         0.00    0.00 48   8.57
 20  8  6  18.51  37.01   0.2429 16.9531   100.75  20  8  6       0.0              0.0         0.00         0.00    0.00 48   8.57
 18 12  6  18.58  37.17   0.2419 17.0888   101.56  18 12  6       0.0              0.1        -0.14         0.00  180.00 48   8.49
 22  4  2  18.58  37.17   0.2419 17.0888   101.56  22  4  2       0.0              0.1        -0.14         0.00  180.00 48   8.49
 20 10  2  18.58  37.17   0.2419 17.0888   101.56  20 10  2       0.0              0.1         0.14         0.00    0.00 48   8.49
 17 13  7  18.64  37.28   0.2412 17.1905   102.16  17 13  7       0.0              0.1        -0.10         0.00  180.00 48   8.44
 19 11  5  18.64  37.28   0.2412 17.1905   102.16  19 11  5       0.0              0.1        -0.10         0.00  180.00 48   8.44
 13 13 13  18.64  37.28   0.2412 17.1905   102.16  13 13 13       0.0              0.1         0.10         0.00    0.00  8   8.44
 16 16  0  18.74  37.47   0.2400 17.3600   103.17  16 16  0       0.0              0.1         0.13         0.00    0.00 12   8.34
 17 15  1  18.79  37.58   0.2393 17.4617   103.77  17 15  1       0.0              0.1         0.09         0.00    0.00 48   8.29
 21  7  5  18.79  37.58   0.2393 17.4617   103.77  21  7  5       0.0              0.1         0.09         0.00    0.00 48   8.29
 15 13 11  18.79  37.58   0.2393 17.4617   103.77  15 13 11       0.0              0.1         0.09         0.00    0.00 48   8.29
 16 14  8  18.81  37.62   0.2391 17.4957   103.97  16 14  8       0.0              0.0         0.00         0.00    0.00 48   8.27
 22  4  4  18.81  37.62   0.2391 17.4957   103.97  22  4  4       0.0              0.0         0.00         0.00    0.00 24   8.27
 20 10  4  18.81  37.62   0.2391 17.4957   103.97  20 10  4       0.0              0.0         0.00         0.00    0.00 48   8.27
 16 16  2  18.81  37.62   0.2391 17.4957   103.97  16 16  2       0.0              0.0         0.00         0.00    0.00 24   8.27
 18 14  0  18.89  37.77   0.2382 17.6313   104.78  18 14  0       0.0              0.1         0.12         0.00    0.00 24   8.19
 22  6  0  18.89  37.77   0.2382 17.6313   104.78  22  6  0       0.0              0.1        -0.12         0.00  180.00 24   8.19
 21  9  1  18.94  37.89   0.2375 17.7330   105.38  21  9  1       0.0              0.1         0.08         0.00    0.00 48   8.14
 19  9  9  18.94  37.89   0.2375 17.7330   105.38  19  9  9       0.0              0.1        -0.08         0.00  180.00 24   8.14
 17 15  3  18.94  37.89   0.2375 17.7330   105.38  17 15  3       0.0              0.1        -0.08         0.00  180.00 48   8.14
 22  6  2  18.96  37.92   0.2372 17.7669   105.59  22  6  2       0.0              0.0         0.00         0.00    0.00 48   8.12
 18 14  2  18.96  37.92   0.2372 17.7669   105.59  18 14  2       0.0              0.0         0.00         0.00    0.00 48   8.12
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 18 10 10  18.96  37.92   0.2372 17.7669   105.59  18 10 10       0.0              0.0         0.00         0.00    0.00 24   8.12
 20  8  8  19.04  38.07   0.2363 17.9025   106.39  20  8  8       0.0              0.1        -0.11         0.00  180.00 24   8.05
 16 16  4  19.04  38.07   0.2363 17.9025   106.39  16 16  4       0.0              0.1        -0.11         0.00  180.00 24   8.05
 17 11 11  19.09  38.19   0.2357 18.0042   107.00  17 11 11       0.0              0.1         0.08         0.00    0.00 24   8.00
 19 11  7  19.09  38.19   0.2357 18.0042   107.00  19 11  7       0.0              0.1        -0.08         0.00  180.00 48   8.00
 23  1  1  19.09  38.19   0.2357 18.0042   107.00  23  1  1       0.0              0.1         0.08         0.00    0.00 24   8.00
 15 15  9  19.09  38.19   0.2357 18.0042   107.00  15 15  9       0.0              0.1         0.08         0.00    0.00 24   8.00
 19 13  1  19.09  38.19   0.2357 18.0042   107.00  19 13  1       0.0              0.1         0.08         0.00    0.00 48   8.00
 21  9  3  19.09  38.19   0.2357 18.0042   107.00  21  9  3       0.0              0.1         0.08         0.00    0.00 48   8.00
 18 12  8  19.11  38.22   0.2355 18.0382   107.20  18 12  8       0.0              0.0         0.00         0.00    0.00 48   7.98
 20 10  6  19.19  38.37   0.2346 18.1738   108.00  20 10  6       0.0              0.1        -0.10         0.00  180.00 48   7.91
 18 14  4  19.19  38.37   0.2346 18.1738   108.00  18 14  4       0.0              0.1        -0.10         0.00  180.00 48   7.91
 14 14 12  19.19  38.37   0.2346 18.1738   108.00  14 14 12       0.0              0.1         0.10         0.00    0.00 24   7.91
 22  6  4  19.19  38.37   0.2346 18.1738   108.00  22  6  4       0.0              0.1         0.10         0.00    0.00 48   7.91
 17 15  5  19.24  38.48   0.2339 18.2755   108.61  17 15  5       0.0              0.1        -0.07         0.00  180.00 48   7.86
 23  3  1  19.24  38.48   0.2339 18.2755   108.61  23  3  1       0.0              0.1        -0.07         0.00  180.00 48   7.86
 21  7  7  19.24  38.48   0.2339 18.2755   108.61  21  7  7       0.0              0.1        -0.07         0.00  180.00 24   7.86
 19 13  3  19.24  38.48   0.2339 18.2755   108.61  19 13  3       0.0              0.1        -0.07         0.00  180.00 48   7.86
 17 13  9  19.24  38.48   0.2339 18.2755   108.61  17 13  9       0.0              0.1         0.07         0.00    0.00 48   7.86
 20 12  0  19.33  38.67   0.2328 18.4450   109.62  20 12  0       0.0              0.1         0.10         0.00    0.00 24   7.78
 16 12 12  19.33  38.67   0.2328 18.4450   109.62  16 12 12       0.0              0.1         0.10         0.00    0.00 24   7.78
 23  3  3  19.39  38.78   0.2322 18.5467   110.22  23  3  3       0.0              0.1        -0.07         0.00  180.00 24   7.73
 21  9  5  19.39  38.78   0.2322 18.5467   110.22  21  9  5       0.0              0.1        -0.07         0.00  180.00 48   7.73
 16 16  6  19.41  38.82   0.2320 18.5807   110.42  16 16  6       0.0              0.0         0.00         0.00    0.00 24   7.72
 20 12  2  19.41  38.82   0.2320 18.5807   110.42  20 12  2       0.0              0.0         0.00         0.00    0.00 48   7.72
 16 14 10  19.48  38.96   0.2311 18.7163   111.23  16 14 10       0.0              0.1         0.09         0.00    0.00 48   7.65
 22  8  2  19.48  38.96   0.2311 18.7163   111.23  22  8  2       0.0              0.1         0.09         0.00    0.00 48   7.65
 23  5  1  19.54  39.07   0.2305 18.8180   111.83  23  5  1       0.0              0.1        -0.06         0.00  180.00 48   7.60
 19 13  5  19.54  39.07   0.2305 18.8180   111.83  19 13  5       0.0              0.1        -0.06         0.00  180.00 48   7.60
 18 14  6  19.56  39.11   0.2303 18.8519   112.03  18 14  6       0.0              0.0         0.00         0.00    0.00 48   7.59
 22  6  6  19.56  39.11   0.2303 18.8519   112.03  22  6  6       0.0              0.0         0.00         0.00    0.00 24   7.59
 20 12  4  19.63  39.26   0.2295 18.9875   112.84  20 12  4       0.0              0.1        -0.08         0.00  180.00 48   7.53
 17 15  7  19.68  39.37   0.2289 19.0892   113.44  17 15  7       0.0              0.1         0.06         0.00    0.00 48   7.48
 23  5  3  19.68  39.37   0.2289 19.0892   113.44  23  5  3       0.0              0.1         0.06         0.00    0.00 48   7.48
 19 11  9  19.68  39.37   0.2289 19.0892   113.44  19 11  9       0.0              0.1         0.06         0.00    0.00 48   7.48
 21 11  1  19.68  39.37   0.2289 19.0892   113.44  21 11  1       0.0              0.1         0.06         0.00    0.00 48   7.48
 15 13 13  19.68  39.37   0.2289 19.0892   113.44  15 13 13       0.0              0.1         0.06         0.00    0.00 24   7.48
 20 10  8  19.70  39.40   0.2287 19.1232   113.65  20 10  8       0.0              0.0         0.00         0.00    0.00 48   7.46
 22  8  4  19.70  39.40   0.2287 19.1232   113.65  22  8  4       0.0              0.0         0.00         0.00    0.00 48   7.46
 18 12 10  19.77  39.55   0.2279 19.2588   114.45  18 12 10       0.0              0.1         0.08         0.00    0.00 48   7.40
 21  9  7  19.83  39.66   0.2273 19.3605   115.06  21  9  7       0.0              0.1        -0.05         0.00  180.00 48   7.36
 21 11  3  19.83  39.66   0.2273 19.3605   115.06  21 11  3       0.0              0.1        -0.05         0.00  180.00 48   7.36
 15 15 11  19.83  39.66   0.2273 19.3605   115.06  15 15 11       0.0              0.1         0.05         0.00    0.00 24   7.36
 16 16  8  19.92  39.84   0.2263 19.5300   116.06  16 16  8       0.0              0.1         0.07         0.00    0.00 24   7.28
 24  0  0  19.92  39.84   0.2263 19.5300   116.06  24  0  0       0.0              0.1         0.07         0.00    0.00  6   7.28
 17 17  1  19.97  39.94   0.2257 19.6317   116.67  17 17  1       0.0              0.1        -0.05         0.00  180.00 24   7.24
 23  5  5  19.97  39.94   0.2257 19.6317   116.67  23  5  5       0.0              0.1         0.05         0.00    0.00 24   7.24
 17 13 11  19.97  39.94   0.2257 19.6317   116.67  17 13 11       0.0              0.1         0.05         0.00    0.00 48   7.24
 23  7  1  19.97  39.94   0.2257 19.6317   116.67  23  7  1       0.0              0.1         0.05         0.00    0.00 48   7.24
 19 13  7  19.97  39.94   0.2257 19.6317   116.67  19 13  7       0.0              0.1         0.05         0.00    0.00 48   7.24
 20 12  6  19.99  39.98   0.2255 19.6657   116.87  20 12  6       0.0              0.0         0.00         0.00    0.00 48   7.23
 18 16  2  20.06  40.12   0.2247 19.8013   117.68  18 16  2       0.0              0.1        -0.07         0.00  180.00 48   7.17
 22 10  0  20.06  40.12   0.2247 19.8013   117.68  22 10  0       0.0              0.1         0.07         0.00    0.00 24   7.17
 22  8  6  20.06  40.12   0.2247 19.8013   117.68  22  8  6       0.0              0.1        -0.07         0.00  180.00 48   7.17
 18 14  8  20.06  40.12   0.2247 19.8013   117.68  18 14  8       0.0              0.1         0.07         0.00    0.00 48   7.17
 24  2  2  20.06  40.12   0.2247 19.8013   117.68  24  2  2       0.0              0.1        -0.07         0.00  180.00 24   7.17
 19 15  1  20.12  40.23   0.2242 19.9030   118.28  19 15  1       0.0              0.0        -0.05         0.00  180.00 48   7.13
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 17 17  3  20.12  40.23   0.2242 19.9030   118.28  17 17  3       0.0              0.0        -0.05         0.00  180.00 24   7.13
 21 11  5  20.12  40.23   0.2242 19.9030   118.28  21 11  5       0.0              0.0        -0.05         0.00  180.00 48   7.13
 23  7  3  20.12  40.23   0.2242 19.9030   118.28  23  7  3       0.0              0.0         0.05         0.00    0.00 48   7.13
 22 10  2  20.13  40.27   0.2240 19.9369   118.48  22 10  2       0.0              0.0         0.00         0.00    0.00 48   7.11
 14 14 14  20.13  40.27   0.2240 19.9369   118.48  14 14 14       0.0              0.0         0.00         0.00    0.00  8   7.11
 24  4  0  20.21  40.41   0.2232 20.0725   119.29  24  4  0       0.0              0.1        -0.06         0.00  180.00 24   7.06
 19 15  3  20.26  40.52   0.2226 20.1742   119.89  19 15  3       0.0              0.0        -0.04         0.00  180.00 48   7.01
 17 15  9  20.26  40.52   0.2226 20.1742   119.89  17 15  9       0.0              0.0         0.04         0.00    0.00 48   7.01
 18 16  4  20.28  40.55   0.2225 20.2082   120.09  18 16  4       0.0              0.0         0.00         0.00    0.00 48   7.00
 24  4  2  20.28  40.55   0.2225 20.2082   120.09  24  4  2       0.0              0.0         0.00         0.00    0.00 48   7.00
 16 14 12  20.28  40.55   0.2225 20.2082   120.09  16 14 12       0.0              0.0         0.00         0.00    0.00 48   7.00
 20 14  2  20.35  40.69   0.2217 20.3438   120.90  20 14  2       0.0              0.1        -0.06         0.00  180.00 48   6.95
 20 10 10  20.35  40.69   0.2217 20.3438   120.90  20 10 10       0.0              0.1         0.06         0.00    0.00 24   6.95
 22 10  4  20.35  40.69   0.2217 20.3438   120.90  22 10  4       0.0              0.1        -0.06         0.00  180.00 48   6.95
 21  9  9  20.40  40.80   0.2212 20.4455   121.50  21  9  9       0.0              0.0         0.04         0.00    0.00 24   6.91
 17 17  5  20.40  40.80   0.2212 20.4455   121.50  17 17  5       0.0              0.0         0.04         0.00    0.00 24   6.91
 23  7  5  20.40  40.80   0.2212 20.4455   121.50  23  7  5       0.0              0.0        -0.04         0.00  180.00 48   6.91
 19 11 11  20.40  40.80   0.2212 20.4455   121.50  19 11 11       0.0              0.0         0.04         0.00    0.00 24   6.91
 24  4  4  20.49  40.98   0.2202 20.6150   122.51  24  4  4       0.0              0.1         0.06         0.00    0.00 24   6.84
 20 12  8  20.49  40.98   0.2202 20.6150   122.51  20 12  8       0.0              0.1         0.06         0.00    0.00 48   6.84
 19 15  5  20.54  41.08   0.2197 20.7167   123.12  19 15  5       0.0              0.0         0.04         0.00    0.00 48   6.80
 23  9  1  20.54  41.08   0.2197 20.7167   123.12  23  9  1       0.0              0.0         0.04         0.00    0.00 48   6.80
 21 13  1  20.54  41.08   0.2197 20.7167   123.12  21 13  1       0.0              0.0        -0.04         0.00  180.00 48   6.80
 19 13  9  20.54  41.08   0.2197 20.7167   123.12  19 13  9       0.0              0.0         0.04         0.00    0.00 48   6.80
 21 11  7  20.54  41.08   0.2197 20.7167   123.12  21 11  7       0.0              0.0         0.04         0.00    0.00 48   6.80
 20 14  4  20.56  41.12   0.2195 20.7507   123.32  20 14  4       0.0              0.0         0.00         0.00    0.00 48   6.79
 16 16 10  20.56  41.12   0.2195 20.7507   123.32  16 16 10       0.0              0.0         0.00         0.00    0.00 24   6.79
 18 12 12  20.56  41.12   0.2195 20.7507   123.32  18 12 12       0.0              0.0         0.00         0.00    0.00 24   6.79
 22  8  8  20.56  41.12   0.2195 20.7507   123.32  22  8  8       0.0              0.0         0.00         0.00    0.00 24   6.79
 18 16  6  20.63  41.26   0.2188 20.8863   124.12  18 16  6       0.0              0.1         0.05         0.00    0.00 48   6.74
 24  6  2  20.63  41.26   0.2188 20.8863   124.12  24  6  2       0.0              0.1         0.05         0.00    0.00 48   6.74
 23  9  3  20.68  41.36   0.2183 20.9880   124.73  23  9  3       0.0              0.0        -0.04         0.00  180.00 48   6.70
 15 15 13  20.68  41.36   0.2183 20.9880   124.73  15 15 13       0.0              0.0        -0.04         0.00  180.00 24   6.70
 21 13  3  20.68  41.36   0.2183 20.9880   124.73  21 13  3       0.0              0.0        -0.04         0.00  180.00 48   6.70
 22 10  6  20.70  41.40   0.2181 21.0219   124.93  22 10  6       0.0              0.0         0.00         0.00    0.00 48   6.69
 18 14 10  20.70  41.40   0.2181 21.0219   124.93  18 14 10       0.0              0.0         0.00         0.00    0.00 48   6.69
 17 17  7  20.82  41.64   0.2169 21.2593   126.34  17 17  7       0.0              0.0         0.03         0.00    0.00 24   6.60
 17 13 13  20.82  41.64   0.2169 21.2593   126.34  17 13 13       0.0              0.0        -0.03         0.00  180.00 24   6.60
 23  7  7  20.82  41.64   0.2169 21.2593   126.34  23  7  7       0.0              0.0        -0.03         0.00  180.00 24   6.60
 25  1  1  20.82  41.64   0.2169 21.2593   126.34  25  1  1       0.0              0.0        -0.03         0.00  180.00 24   6.60
 24  6  4  20.84  41.68   0.2167 21.2932   126.54  24  6  4       0.0              0.0         0.00         0.00    0.00 48   6.59
 22 12  2  20.91  41.81   0.2160 21.4288   127.35  22 12  2       0.0              0.0        -0.04         0.00  180.00 48   6.54
 20 14  6  20.91  41.81   0.2160 21.4288   127.35  20 14  6       0.0              0.0         0.04         0.00    0.00 48   6.54
 25  3  1  20.96  41.92   0.2155 21.5305   127.95  25  3  1       0.0              0.0        -0.03         0.00  180.00 48   6.50
 19 15  7  20.96  41.92   0.2155 21.5305   127.95  19 15  7       0.0              0.0         0.03         0.00    0.00 48   6.50
 23  9  5  20.96  41.92   0.2155 21.5305   127.95  23  9  5       0.0              0.0        -0.03         0.00  180.00 48   6.50
 17 15 11  20.96  41.92   0.2155 21.5305   127.95  17 15 11       0.0              0.0        -0.03         0.00  180.00 48   6.50
 21 13  5  20.96  41.92   0.2155 21.5305   127.95  21 13  5       0.0              0.0         0.03         0.00    0.00 48   6.50
 24  8  0  21.05  42.09   0.2147 21.7000   128.96  24  8  0       0.0              0.0         0.04         0.00    0.00 24   6.44
 21 11  9  21.10  42.19   0.2142 21.8018   129.56  21 11  9       0.0              0.0         0.03         0.00    0.00 48   6.41
 25  3  3  21.10  42.19   0.2142 21.8018   129.56  25  3  3       0.0              0.0         0.03         0.00    0.00 24   6.41
 24  8  2  21.11  42.23   0.2140 21.8357   129.77  24  8  2       0.0              0.0         0.00         0.00    0.00 48   6.40
 22 12  4  21.11  42.23   0.2140 21.8357   129.77  22 12  4       0.0              0.0         0.00         0.00    0.00 48   6.40
 20 12 10  21.11  42.23   0.2140 21.8357   129.77  20 12 10       0.0              0.0         0.00         0.00    0.00 48   6.40
 18 16  8  21.11  42.23   0.2140 21.8357   129.77  18 16  8       0.0              0.0         0.00         0.00    0.00 48   6.40
 22 10  8  21.18  42.37   0.2133 21.9713   130.57  22 10  8       0.0              0.0         0.04         0.00    0.00 48   6.35
 18 18  0  21.18  42.37   0.2133 21.9713   130.57  18 18  0       0.0              0.0        -0.04         0.00  180.00 12   6.35
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 24  6  6  21.18  42.37   0.2133 21.9713   130.57  24  6  6       0.0              0.0        -0.04         0.00  180.00 24   6.35
 16 14 14  21.18  42.37   0.2133 21.9713   130.57  16 14 14       0.0              0.0        -0.04         0.00  180.00 24   6.35
 25  5  1  21.23  42.47   0.2128 22.0730   131.18  25  5  1       0.0              0.0         0.03         0.00    0.00 48   6.31
 19 13 11  21.23  42.47   0.2128 22.0730   131.18  19 13 11       0.0              0.0        -0.03         0.00  180.00 48   6.31
 19 17  1  21.23  42.47   0.2128 22.0730   131.18  19 17  1       0.0              0.0        -0.03         0.00  180.00 48   6.31
 23 11  1  21.23  42.47   0.2128 22.0730   131.18  23 11  1       0.0              0.0        -0.03         0.00  180.00 48   6.31
 18 18  2  21.25  42.50   0.2127 22.1069   131.38  18 18  2       0.0              0.0         0.00         0.00    0.00 24   6.30
 24  8  4  21.32  42.64   0.2120 22.2425   132.18  24  8  4       0.0              0.0        -0.04         0.00  180.00 48   6.26
 20 16  0  21.32  42.64   0.2120 22.2425   132.18  20 16  0       0.0              0.0        -0.04         0.00  180.00 24   6.26
 16 16 12  21.32  42.64   0.2120 22.2425   132.18  16 16 12       0.0              0.0        -0.04         0.00  180.00 24   6.26
 23  9  7  21.37  42.74   0.2116 22.3443   132.79  23  9  7       0.0              0.0         0.03         0.00    0.00 48   6.22
 19 17  3  21.37  42.74   0.2116 22.3443   132.79  19 17  3       0.0              0.0         0.03         0.00    0.00 48   6.22
 17 17  9  21.37  42.74   0.2116 22.3443   132.79  17 17  9       0.0              0.0        -0.03         0.00  180.00 24   6.22
 23 11  3  21.37  42.74   0.2116 22.3443   132.79  23 11  3       0.0              0.0        -0.03         0.00  180.00 48   6.22
 21 13  7  21.37  42.74   0.2116 22.3443   132.79  21 13  7       0.0              0.0         0.03         0.00    0.00 48   6.22
 25  5  3  21.37  42.74   0.2116 22.3443   132.79  25  5  3       0.0              0.0         0.03         0.00    0.00 48   6.22
 20 14  8  21.39  42.78   0.2114 22.3782   132.99  20 14  8       0.0              0.0         0.00         0.00    0.00 48   6.21
 20 16  2  21.39  42.78   0.2114 22.3782   132.99  20 16  2       0.0              0.0         0.00         0.00    0.00 48   6.21
 22 12  6  21.46  42.91   0.2108 22.5138   133.80  22 12  6       0.0              0.0         0.03         0.00    0.00 48   6.17
 18 14 12  21.46  42.91   0.2108 22.5138   133.80  18 14 12       0.0              0.0        -0.03         0.00  180.00 48   6.17
 18 18  4  21.46  42.91   0.2108 22.5138   133.80  18 18  4       0.0              0.0         0.03         0.00    0.00 24   6.17
 19 15  9  21.51  43.01   0.2103 22.6155   134.40  19 15  9       0.0              0.0        -0.02         0.00  180.00 48   6.14
 21 15  1  21.51  43.01   0.2103 22.6155   134.40  21 15  1       0.0              0.0        -0.02         0.00  180.00 48   6.14
 20 16  4  21.59  43.18   0.2095 22.7850   135.41  20 16  4       0.0              0.0         0.03         0.00    0.00 48   6.08
 25  7  1  21.64  43.28   0.2090 22.8868   136.01  25  7  1       0.0              0.0         0.02         0.00    0.00 48   6.05
 23 11  5  21.64  43.28   0.2090 22.8868   136.01  23 11  5       0.0              0.0         0.02         0.00    0.00 48   6.05
 25  5  5  21.64  43.28   0.2090 22.8868   136.01  25  5  5       0.0              0.0        -0.02         0.00  180.00 24   6.05
 19 17  5  21.64  43.28   0.2090 22.8868   136.01  19 17  5       0.0              0.0         0.02         0.00    0.00 48   6.05
 21 15  3  21.64  43.28   0.2090 22.8868   136.01  21 15  3       0.0              0.0         0.02         0.00    0.00 48   6.05
 15 15 15  21.64  43.28   0.2090 22.8868   136.01  15 15 15       0.0              0.0        -0.02         0.00  180.00  8   6.05
 24  8  6  21.66  43.32   0.2089 22.9207   136.21  24  8  6       0.0              0.0         0.00         0.00    0.00 48   6.04
 24 10  2  21.73  43.45   0.2083 23.0563   137.02  24 10  2       0.0              0.0        -0.03         0.00  180.00 48   6.00
 18 16 10  21.73  43.45   0.2083 23.0563   137.02  18 16 10       0.0              0.0        -0.03         0.00  180.00 48   6.00
 26  2  0  21.73  43.45   0.2083 23.0563   137.02  26  2  0       0.0              0.0        -0.03         0.00  180.00 24   6.00
 22 14  0  21.73  43.45   0.2083 23.0563   137.02  22 14  0       0.0              0.0        -0.03         0.00  180.00 24   6.00
 25  7  3  21.78  43.55   0.2078 23.1580   137.62  25  7  3       0.0              0.0        -0.02         0.00  180.00 48   5.97
 21 11 11  21.78  43.55   0.2078 23.1580   137.62  21 11 11       0.0              0.0        -0.02         0.00  180.00 24   5.97
 17 15 13  21.78  43.55   0.2078 23.1580   137.62  17 15 13       0.0              0.0        -0.02         0.00  180.00 48   5.97
 22 14  2  21.79  43.59   0.2076 23.1919   137.83  22 14  2       0.0              0.0         0.00         0.00    0.00 48   5.96
 22 10 10  21.79  43.59   0.2076 23.1919   137.83  22 10 10       0.0              0.0         0.00         0.00    0.00 24   5.96
 26  2  2  21.79  43.59   0.2076 23.1919   137.83  26  2  2       0.0              0.0         0.00         0.00    0.00 24   5.96
 18 18  6  21.79  43.59   0.2076 23.1919   137.83  18 18  6       0.0              0.0         0.00         0.00    0.00 24   5.96
 20 12 12  21.86  43.72   0.2070 23.3275   138.63  20 12 12       0.0              0.0        -0.03         0.00  180.00 24   5.92
 21 13  9  21.91  43.82   0.2066 23.4293   139.24  21 13  9       0.0              0.0        -0.02         0.00  180.00 48   5.89
 21 15  5  21.91  43.82   0.2066 23.4293   139.24  21 15  5       0.0              0.0         0.02         0.00    0.00 48   5.89
 23  9  9  21.91  43.82   0.2066 23.4293   139.24  23  9  9       0.0              0.0         0.02         0.00    0.00 24   5.89
 20 16  6  21.93  43.85   0.2064 23.4632   139.44  20 16  6       0.0              0.0         0.00         0.00    0.00 48   5.88
 24 10  4  21.93  43.85   0.2064 23.4632   139.44  24 10  4       0.0              0.0         0.00         0.00    0.00 48   5.88
 22 12  8  21.93  43.85   0.2064 23.4632   139.44  22 12  8       0.0              0.0         0.00         0.00    0.00 48   5.88
 26  4  2  21.99  43.99   0.2059 23.5988   140.24  26  4  2       0.0              0.0         0.03         0.00    0.00 48   5.84
 20 14 10  21.99  43.99   0.2059 23.5988   140.24  20 14 10       0.0              0.0        -0.03         0.00  180.00 48   5.84
 22 14  4  21.99  43.99   0.2059 23.5988   140.24  22 14  4       0.0              0.0         0.03         0.00    0.00 48   5.84
 19 13 13  22.04  44.09   0.2054 23.7005   140.85  19 13 13       0.0              0.0        -0.02         0.00  180.00 24   5.81
 23 11  7  22.04  44.09   0.2054 23.7005   140.85  23 11  7       0.0              0.0         0.02         0.00    0.00 48   5.81
 19 17  7  22.04  44.09   0.2054 23.7005   140.85  19 17  7       0.0              0.0        -0.02         0.00  180.00 48   5.81
 17 17 11  22.04  44.09   0.2054 23.7005   140.85  17 17 11       0.0              0.0        -0.02         0.00  180.00 24   5.81
 25  7  5  22.04  44.09   0.2054 23.7005   140.85  25  7  5       0.0              0.0        -0.02         0.00  180.00 48   5.81
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 23 13  1  22.04  44.09   0.2054 23.7005   140.85  23 13  1       0.0              0.0        -0.02         0.00  180.00 48   5.81
 24  8  8  22.13  44.25   0.2047 23.8700   141.86  24  8  8       0.0              0.0         0.02         0.00    0.00 24   5.76
 25  9  1  22.18  44.35   0.2042 23.9718   142.46  25  9  1       0.0              0.0        -0.02         0.00  180.00 48   5.73
 19 15 11  22.18  44.35   0.2042 23.9718   142.46  19 15 11       0.0              0.0        -0.02         0.00  180.00 48   5.73
 23 13  3  22.18  44.35   0.2042 23.9718   142.46  23 13  3       0.0              0.0         0.02         0.00    0.00 48   5.73
 16 16 14  22.19  44.38   0.2041 24.0057   142.66  16 16 14       0.0              0.0         0.00         0.00    0.00 24   5.72
 26  4  4  22.19  44.38   0.2041 24.0057   142.66  26  4  4       0.0              0.0         0.00         0.00    0.00 24   5.72
 18 18  8  22.26  44.52   0.2035 24.1413   143.47  18 18  8       0.0              0.0        -0.02         0.00  180.00 24   5.68
 26  6  0  22.26  44.52   0.2035 24.1413   143.47  26  6  0       0.0              0.0         0.02         0.00    0.00 24   5.68
 24 10  6  22.26  44.52   0.2035 24.1413   143.47  24 10  6       0.0              0.0         0.02         0.00    0.00 48   5.68
 21 15  7  22.31  44.61   0.2031 24.2430   144.07  21 15  7       0.0              0.0        -0.02         0.00  180.00 48   5.65
 25  9  3  22.31  44.61   0.2031 24.2430   144.07  25  9  3       0.0              0.0        -0.02         0.00  180.00 48   5.65
 18 14 14  22.32  44.65   0.2030 24.2769   144.27  18 14 14       0.0              0.0         0.00         0.00    0.00 24   5.64
 26  6  2  22.32  44.65   0.2030 24.2769   144.27  26  6  2       0.0              0.0         0.00         0.00    0.00 48   5.64
 22 14  6  22.32  44.65   0.2030 24.2769   144.27  22 14  6       0.0              0.0         0.00         0.00    0.00 48   5.64
 24 12  0  22.39  44.78   0.2024 24.4125   145.08  24 12  0       0.0              0.0        -0.02         0.00  180.00 24   5.61
 20 16  8  22.39  44.78   0.2024 24.4125   145.08  20 16  8       0.0              0.0        -0.02         0.00  180.00 48   5.61
 25  7  7  22.44  44.88   0.2020 24.5143   145.68  25  7  7       0.0              0.0         0.01         0.00    0.00 24   5.58
 23 13  5  22.44  44.88   0.2020 24.5143   145.68  23 13  5       0.0              0.0         0.01         0.00    0.00 48   5.58
 19 19  1  22.44  44.88   0.2020 24.5143   145.68  19 19  1       0.0              0.0         0.01         0.00    0.00 24   5.58
 18 16 12  22.45  44.91   0.2018 24.5482   145.89  18 16 12       0.0              0.0         0.00         0.00    0.00 48   5.57
 24 12  2  22.45  44.91   0.2018 24.5482   145.89  24 12  2       0.0              0.0         0.00         0.00    0.00 48   5.57
 22 12 10  22.52  45.04   0.2013 24.6838   146.69  22 12 10       0.0              0.0        -0.02         0.00  180.00 48   5.53
 26  6  4  22.52  45.04   0.2013 24.6838   146.69  26  6  4       0.0              0.0        -0.02         0.00  180.00 48   5.53
 20 18  2  22.52  45.04   0.2013 24.6838   146.69  20 18  2       0.0              0.0         0.02         0.00    0.00 48   5.53
 25  9  5  22.57  45.14   0.2009 24.7855   147.30  25  9  5       0.0              0.0         0.01         0.00    0.00 48   5.51
 27  1  1  22.57  45.14   0.2009 24.7855   147.30  27  1  1       0.0              0.0        -0.01         0.00  180.00 24   5.51
 21 13 11  22.57  45.14   0.2009 24.7855   147.30  21 13 11       0.0              0.0        -0.01         0.00  180.00 48   5.51
 19 17  9  22.57  45.14   0.2009 24.7855   147.30  19 17  9       0.0              0.0        -0.01         0.00  180.00 48   5.51
 23 11  9  22.57  45.14   0.2009 24.7855   147.30  23 11  9       0.0              0.0        -0.01         0.00  180.00 48   5.51
 21 17  1  22.57  45.14   0.2009 24.7855   147.30  21 17  1       0.0              0.0         0.01         0.00    0.00 48   5.51
 19 19  3  22.57  45.14   0.2009 24.7855   147.30  19 19  3       0.0              0.0         0.01         0.00    0.00 24   5.51
 24 12  4  22.65  45.30   0.2002 24.9550   148.30  24 12  4       0.0              0.0         0.02         0.00    0.00 48   5.46
 27  3  1  22.70  45.40   0.1998 25.0568   148.91  27  3  1       0.0              0.0         0.01         0.00    0.00 48   5.43
 21 17  3  22.70  45.40   0.1998 25.0568   148.91  21 17  3       0.0              0.0         0.01         0.00    0.00 48   5.43
 17 15 15  22.70  45.40   0.1998 25.0568   148.91  17 15 15       0.0              0.0         0.01         0.00    0.00 24   5.43
 20 14 12  22.72  45.43   0.1996 25.0907   149.11  20 14 12       0.0              0.0         0.00         0.00    0.00 48   5.43
 24 10  8  22.72  45.43   0.1996 25.0907   149.11  24 10  8       0.0              0.0         0.00         0.00    0.00 48   5.43
 20 18  4  22.72  45.43   0.1996 25.0907   149.11  20 18  4       0.0              0.0         0.00         0.00    0.00 48   5.43
 22 14  8  22.78  45.56   0.1991 25.2263   149.92  22 14  8       0.0              0.0        -0.02         0.00  180.00 48   5.39
 26  8  2  22.78  45.56   0.1991 25.2263   149.92  26  8  2       0.0              0.0        -0.02         0.00  180.00 48   5.39
 22 16  2  22.78  45.56   0.1991 25.2263   149.92  22 16  2       0.0              0.0         0.02         0.00    0.00 48   5.39
 27  3  3  22.83  45.66   0.1987 25.3280   150.52  27  3  3       0.0              0.0         0.01         0.00    0.00 24   5.36
 21 15  9  22.83  45.66   0.1987 25.3280   150.52  21 15  9       0.0              0.0        -0.01         0.00  180.00 48   5.36
 25 11  1  22.83  45.66   0.1987 25.3280   150.52  25 11  1       0.0              0.0        -0.01         0.00  180.00 48   5.36
 19 19  5  22.83  45.66   0.1987 25.3280   150.52  19 19  5       0.0              0.0        -0.01         0.00  180.00 24   5.36
 23 13  7  22.83  45.66   0.1987 25.3280   150.52  23 13  7       0.0              0.0        -0.01         0.00  180.00 48   5.36
 17 17 13  22.83  45.66   0.1987 25.3280   150.52  17 17 13       0.0              0.0         0.01         0.00    0.00 24   5.36
 18 18 10  22.84  45.69   0.1986 25.3619   150.72  18 18 10       0.0              0.0         0.00         0.00    0.00 24   5.36
 26  6  6  22.84  45.69   0.1986 25.3619   150.72  26  6  6       0.0              0.0         0.00         0.00    0.00 24   5.36
 25 11  3  22.96  45.91   0.1976 25.5993   152.13  25 11  3       0.0              0.0         0.01         0.00    0.00 48   5.30
 21 17  5  22.96  45.91   0.1976 25.5993   152.13  21 17  5       0.0              0.0        -0.01         0.00  180.00 48   5.30
 27  5  1  22.96  45.91   0.1976 25.5993   152.13  27  5  1       0.0              0.0         0.01         0.00    0.00 48   5.30
 23 15  1  22.96  45.91   0.1976 25.5993   152.13  23 15  1       0.0              0.0         0.01         0.00    0.00 48   5.30
 19 15 13  22.96  45.91   0.1976 25.5993   152.13  19 15 13       0.0              0.0         0.01         0.00    0.00 48   5.30
 25  9  7  22.96  45.91   0.1976 25.5993   152.13  25  9  7       0.0              0.0         0.01         0.00    0.00 48   5.30
 24 12  6  22.97  45.95   0.1975 25.6332   152.33  24 12  6       0.0              0.0         0.00         0.00    0.00 48   5.29
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 22 16  4  22.97  45.95   0.1975 25.6332   152.33  22 16  4       0.0              0.0         0.00         0.00    0.00 48   5.29
 20 16 10  22.97  45.95   0.1975 25.6332   152.33  20 16 10       0.0              0.0         0.00         0.00    0.00 48   5.29
 26  8  4  22.97  45.95   0.1975 25.6332   152.33  26  8  4       0.0              0.0         0.00         0.00    0.00 48   5.29
 20 18  6  23.04  46.07   0.1970 25.7688   153.14  20 18  6       0.0              0.0        -0.01         0.00  180.00 48   5.26
 27  5  3  23.09  46.17   0.1966 25.8705   153.75  27  5  3       0.0              0.0        -0.01         0.00  180.00 48   5.23
 23 15  3  23.09  46.17   0.1966 25.8705   153.75  23 15  3       0.0              0.0         0.01         0.00    0.00 48   5.23
 16 16 16  23.17  46.33   0.1960 26.0400   154.75  16 16 16       0.0              0.0         0.01         0.00    0.00  8   5.19
 23 11 11  23.21  46.43   0.1956 26.1418   155.36  23 11 11       0.0              0.0        -0.01         0.00  180.00 24   5.17
 19 17 11  23.21  46.43   0.1956 26.1418   155.36  19 17 11       0.0              0.0         0.01         0.00    0.00 48   5.17
 19 19  7  23.21  46.43   0.1956 26.1418   155.36  19 19  7       0.0              0.0        -0.01         0.00  180.00 24   5.17
 25 11  5  23.21  46.43   0.1956 26.1418   155.36  25 11  5       0.0              0.0         0.01         0.00    0.00 48   5.17
 22 12 12  23.23  46.46   0.1955 26.1757   155.56  22 12 12       0.0              0.0         0.00         0.00    0.00 24   5.16
 26 10  0  23.29  46.59   0.1950 26.3113   156.36  26 10  0       0.0              0.0        -0.01         0.00  180.00 24   5.13
 24 10 10  23.29  46.59   0.1950 26.3113   156.36  24 10 10       0.0              0.0        -0.01         0.00  180.00 24   5.13
 18 16 14  23.29  46.59   0.1950 26.3113   156.36  18 16 14       0.0              0.0         0.01         0.00    0.00 48   5.13
 24 14  2  23.29  46.59   0.1950 26.3113   156.36  24 14  2       0.0              0.0         0.01         0.00    0.00 48   5.13
 22 16  6  23.29  46.59   0.1950 26.3113   156.36  22 16  6       0.0              0.0        -0.01         0.00  180.00 48   5.13
 26  8  6  23.29  46.59   0.1950 26.3113   156.36  26  8  6       0.0              0.0         0.01         0.00    0.00 48   5.13
 27  7  1  23.34  46.68   0.1946 26.4130   156.97  27  7  1       0.0              0.0        -0.01         0.00  180.00 48   5.10
 21 17  7  23.34  46.68   0.1946 26.4130   156.97  21 17  7       0.0              0.0        -0.01         0.00  180.00 48   5.10
 23 15  5  23.34  46.68   0.1946 26.4130   156.97  23 15  5       0.0              0.0        -0.01         0.00  180.00 48   5.10
 23 13  9  23.34  46.68   0.1946 26.4130   156.97  23 13  9       0.0              0.0        -0.01         0.00  180.00 48   5.10
 27  5  5  23.34  46.68   0.1946 26.4130   156.97  27  5  5       0.0              0.0        -0.01         0.00  180.00 24   5.10
 21 13 13  23.34  46.68   0.1946 26.4130   156.97  21 13 13       0.0              0.0         0.01         0.00    0.00 24   5.10
 26 10  2  23.36  46.71   0.1945 26.4469   157.17  26 10  2       0.0              0.0         0.00         0.00    0.00 48   5.09
 22 14 10  23.36  46.71   0.1945 26.4469   157.17  22 14 10       0.0              0.0         0.00         0.00    0.00 48   5.09
 24 12  8  23.42  46.84   0.1940 26.5825   157.98  24 12  8       0.0              0.0        -0.01         0.00  180.00 48   5.06
 28  0  0  23.42  46.84   0.1940 26.5825   157.98  28  0  0       0.0              0.0        -0.01         0.00  180.00  6   5.06
 27  7  3  23.47  46.93   0.1936 26.6843   158.58  27  7  3       0.0              0.0        -0.01         0.00  180.00 48   5.04
 21 15 11  23.47  46.93   0.1936 26.6843   158.58  21 15 11       0.0              0.0         0.01         0.00    0.00 48   5.04
 25  9  9  23.47  46.93   0.1936 26.6843   158.58  25  9  9       0.0              0.0        -0.01         0.00  180.00 24   5.04
 24 14  4  23.48  46.97   0.1935 26.7182   158.78  24 14  4       0.0              0.0         0.00         0.00    0.00 48   5.03
 20 18  8  23.48  46.97   0.1935 26.7182   158.78  20 18  8       0.0              0.0         0.00         0.00    0.00 48   5.03
 28  2  2  23.55  47.09   0.1930 26.8538   159.59  28  2  2       0.0              0.0         0.01         0.00    0.00 24   5.00
 18 18 12  23.55  47.09   0.1930 26.8538   159.59  18 18 12       0.0              0.0         0.01         0.00    0.00 24   5.00
 20 14 14  23.55  47.09   0.1930 26.8538   159.59  20 14 14       0.0              0.0         0.01         0.00    0.00 24   5.00
 26 10  4  23.55  47.09   0.1930 26.8538   159.59  26 10  4       0.0              0.0         0.01         0.00    0.00 48   5.00
 25 13  1  23.59  47.19   0.1926 26.9555   160.19  25 13  1       0.0              0.0         0.01         0.00    0.00 48   4.98
 25 11  7  23.59  47.19   0.1926 26.9555   160.19  25 11  7       0.0              0.0        -0.01         0.00  180.00 48   4.98
 20 20  0  23.67  47.34   0.1920 27.1250   161.20  20 20  0       0.0              0.0         0.01         0.00    0.00 12   4.94
 20 16 12  23.67  47.34   0.1920 27.1250   161.20  20 16 12       0.0              0.0         0.01         0.00    0.00 48   4.94
 28  4  0  23.67  47.34   0.1920 27.1250   161.20  28  4  0       0.0              0.0         0.01         0.00    0.00 24   4.94
 21 19  1  23.72  47.44   0.1916 27.2268   161.81  21 19  1       0.0              0.0         0.01         0.00    0.00 48   4.92
 19 19  9  23.72  47.44   0.1916 27.2268   161.81  19 19  9       0.0              0.0         0.01         0.00    0.00 24   4.92
 25 13  3  23.72  47.44   0.1916 27.2268   161.81  25 13  3       0.0              0.0         0.01         0.00    0.00 48   4.92
 23 15  7  23.72  47.44   0.1916 27.2268   161.81  23 15  7       0.0              0.0        -0.01         0.00  180.00 48   4.92
 17 17 15  23.72  47.44   0.1916 27.2268   161.81  17 17 15       0.0              0.0         0.01         0.00    0.00 24   4.92
 27  7  5  23.72  47.44   0.1916 27.2268   161.81  27  7  5       0.0              0.0         0.01         0.00    0.00 48   4.92
 20 20  2  23.73  47.47   0.1915 27.2607   162.01  20 20  2       0.0              0.0         0.00         0.00    0.00 24   4.91
 28  4  2  23.73  47.47   0.1915 27.2607   162.01  28  4  2       0.0              0.0         0.00         0.00    0.00 48   4.91
 22 16  8  23.73  47.47   0.1915 27.2607   162.01  22 16  8       0.0              0.0         0.00         0.00    0.00 48   4.91
 26  8  8  23.73  47.47   0.1915 27.2607   162.01  26  8  8       0.0              0.0         0.00         0.00    0.00 24   4.91
 24 14  6  23.80  47.59   0.1911 27.3963   162.81  24 14  6       0.0              0.0        -0.01         0.00  180.00 48   4.88
 22 18  0  23.80  47.59   0.1911 27.3963   162.81  22 18  0       0.0              0.0         0.01         0.00    0.00 24   4.88
 27  9  1  23.84  47.69   0.1907 27.4980   163.42  27  9  1       0.0              0.0        -0.01         0.00  180.00 48   4.86
 19 15 15  23.84  47.69   0.1907 27.4980   163.42  19 15 15       0.0              0.0         0.01         0.00    0.00 24   4.86
 21 19  3  23.84  47.69   0.1907 27.4980   163.42  21 19  3       0.0              0.0        -0.01         0.00  180.00 48   4.86
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 21 17  9  23.84  47.69   0.1907 27.4980   163.42  21 17  9       0.0              0.0         0.01         0.00    0.00 48   4.86
 22 18  2  23.86  47.72   0.1906 27.5319   163.62  22 18  2       0.0              0.0         0.00         0.00    0.00 48   4.85
 26 10  6  23.86  47.72   0.1906 27.5319   163.62  26 10  6       0.0              0.0         0.00         0.00    0.00 48   4.85
 20 20  4  23.92  47.84   0.1901 27.6675   164.42  20 20  4       0.0              0.0        -0.01         0.00  180.00 24   4.83
 28  4  4  23.92  47.84   0.1901 27.6675   164.42  28  4  4       0.0              0.0        -0.01         0.00  180.00 24   4.83
 25 13  5  23.97  47.94   0.1898 27.7693   165.03  25 13  5       0.0              0.0        -0.01         0.00  180.00 48   4.80
 23 13 11  23.97  47.94   0.1898 27.7693   165.03  23 13 11       0.0              0.0         0.01         0.00    0.00 48   4.80
 23 17  1  23.97  47.94   0.1898 27.7693   165.03  23 17  1       0.0              0.0         0.01         0.00    0.00 48   4.80
 19 17 13  23.97  47.94   0.1898 27.7693   165.03  19 17 13       0.0              0.0         0.01         0.00    0.00 48   4.80
 27  9  3  23.97  47.94   0.1898 27.7693   165.03  27  9  3       0.0              0.0         0.01         0.00    0.00 48   4.80
 24 12 10  23.98  47.97   0.1897 27.8032   165.23  24 12 10       0.0              0.0         0.00         0.00    0.00 48   4.80
 22 18  4  24.05  48.09   0.1892 27.9388   166.04  22 18  4       0.0              0.0        -0.01         0.00  180.00 48   4.77
 26 12  2  24.05  48.09   0.1892 27.9388   166.04  26 12  2       0.0              0.0         0.01         0.00    0.00 48   4.77
 20 18 10  24.05  48.09   0.1892 27.9388   166.04  20 18 10       0.0              0.0         0.01         0.00    0.00 48   4.77
 28  6  2  24.05  48.09   0.1892 27.9388   166.04  28  6  2       0.0              0.0        -0.01         0.00  180.00 48   4.77
 22 14 12  24.05  48.09   0.1892 27.9388   166.04  22 14 12       0.0              0.0         0.01         0.00    0.00 48   4.77
 27  7  7  24.09  48.19   0.1888 28.0405   166.64  27  7  7       0.0              0.0         0.01         0.00    0.00 24   4.75
 25 11  9  24.09  48.19   0.1888 28.0405   166.64  25 11  9       0.0              0.0        -0.01         0.00  180.00 48   4.75
 21 19  5  24.09  48.19   0.1888 28.0405   166.64  21 19  5       0.0              0.0        -0.01         0.00  180.00 48   4.75
 23 17  3  24.09  48.19   0.1888 28.0405   166.64  23 17  3       0.0              0.0        -0.01         0.00  180.00 48   4.75
 24 16  0  24.17  48.34   0.1883 28.2100   167.65  24 16  0       0.0              0.0         0.01         0.00    0.00 24   4.71
 27  9  5  24.22  48.43   0.1879 28.3118   168.25  27  9  5       0.0              0.0         0.01         0.00    0.00 48   4.69
 21 15 13  24.22  48.43   0.1879 28.3118   168.25  21 15 13       0.0              0.0         0.01         0.00    0.00 48   4.69
 23 15  9  24.22  48.43   0.1879 28.3118   168.25  23 15  9       0.0              0.0         0.01         0.00    0.00 48   4.69
 24 14  8  24.23  48.46   0.1878 28.3457   168.45  24 14  8       0.0              0.0         0.00         0.00    0.00 48   4.69
 26 12  4  24.23  48.46   0.1878 28.3457   168.45  26 12  4       0.0              0.0         0.00         0.00    0.00 48   4.69
 24 16  2  24.23  48.46   0.1878 28.3457   168.45  24 16  2       0.0              0.0         0.00         0.00    0.00 48   4.69
 20 20  6  24.23  48.46   0.1878 28.3457   168.45  20 20  6       0.0              0.0         0.00         0.00    0.00 24   4.69
 28  6  4  24.23  48.46   0.1878 28.3457   168.45  28  6  4       0.0              0.0         0.00         0.00    0.00 48   4.69
 18 16 16  24.23  48.46   0.1878 28.3457   168.45  18 16 16       0.0              0.0         0.00         0.00    0.00 24   4.69
 22 16 10  24.29  48.59   0.1874 28.4813   169.26  22 16 10       0.0              0.0         0.01         0.00    0.00 48   4.66
 26 10  8  24.29  48.59   0.1874 28.4813   169.26  26 10  8       0.0              0.0        -0.01         0.00  180.00 48   4.66
 29  1  1  24.34  48.68   0.1870 28.5830   169.87  29  1  1       0.0              0.0         0.01         0.00    0.00 24   4.64
 19 19 11  24.34  48.68   0.1870 28.5830   169.87  19 19 11       0.0              0.0         0.01         0.00    0.00 24   4.64
 23 17  5  24.34  48.68   0.1870 28.5830   169.87  23 17  5       0.0              0.0        -0.01         0.00  180.00 48   4.64
 25 13  7  24.34  48.68   0.1870 28.5830   169.87  25 13  7       0.0              0.0        -0.01         0.00  180.00 48   4.64
 22 18  6  24.36  48.71   0.1869 28.6169   170.07  22 18  6       0.0              0.0         0.00         0.00    0.00 48   4.63
 18 18 14  24.36  48.71   0.1869 28.6169   170.07  18 18 14       0.0              0.0         0.00         0.00    0.00 24   4.63
 24 16  4  24.42  48.83   0.1865 28.7525   170.87  24 16  4       0.0              0.0        -0.01         0.00  180.00 48   4.61
 28  8  0  24.42  48.83   0.1865 28.7525   170.87  28  8  0       0.0              0.0        -0.01         0.00  180.00 24   4.61
 21 19  7  24.46  48.93   0.1862 28.8543   171.48  21 19  7       0.0              0.0         0.00         0.00    0.00 48   4.59
 25 15  1  24.46  48.93   0.1862 28.8543   171.48  25 15  1       0.0              0.0         0.00         0.00    0.00 48   4.59
 29  3  1  24.46  48.93   0.1862 28.8543   171.48  29  3  1       0.0              0.0         0.00         0.00    0.00 48   4.59
 27 11  1  24.46  48.93   0.1862 28.8543   171.48  27 11  1       0.0              0.0         0.00         0.00    0.00 48   4.59
 21 17 11  24.46  48.93   0.1862 28.8543   171.48  21 17 11       0.0              0.0         0.00         0.00    0.00 48   4.59
 20 16 14  24.48  48.96   0.1861 28.8882   171.68  20 16 14       0.0              0.0         0.00         0.00    0.00 48   4.58
 28  8  2  24.48  48.96   0.1861 28.8882   171.68  28  8  2       0.0              0.0         0.00         0.00    0.00 48   4.58
 28  6  6  24.54  49.08   0.1856 29.0238   172.48  28  6  6       0.0              0.0         0.01         0.00    0.00 24   4.55
 26 12  6  24.54  49.08   0.1856 29.0238   172.48  26 12  6       0.0              0.0        -0.01         0.00  180.00 48   4.55
 29  3  3  24.58  49.17   0.1853 29.1255   173.09  29  3  3       0.0              0.0         0.00         0.00  180.00 24   4.53
 27  9  7  24.58  49.17   0.1853 29.1255   173.09  27  9  7       0.0              0.0         0.00         0.00  180.00 48   4.53
 27 11  3  24.58  49.17   0.1853 29.1255   173.09  27 11  3       0.0              0.0         0.00         0.00    0.00 48   4.53
 25 15  3  24.58  49.17   0.1853 29.1255   173.09  25 15  3       0.0              0.0         0.00         0.00  180.00 48   4.53
 28  8  4  24.66  49.32   0.1848 29.2950   174.10  28  8  4       0.0              0.0         0.01         0.00    0.00 48   4.50
 24 12 12  24.66  49.32   0.1848 29.2950   174.10  24 12 12       0.0              0.0         0.01         0.00    0.00 24   4.50
 20 20  8  24.66  49.32   0.1848 29.2950   174.10  20 20  8       0.0              0.0         0.01         0.00    0.00 24   4.50
 25 11 11  24.71  49.41   0.1844 29.3968   174.70  25 11 11       0.0              0.0         0.00         0.00    0.00 24   4.48
1

  H  K  L  THETA  2THETA D VALUE  1/D**2 SIN2*1000  H  K  L INTENSITY         /F(HKL)/       A(HKL)      B(HKL) PHA.ANG. MULT   LPG

 29  5  1  24.71  49.41   0.1844 29.3968   174.70  29  5  1       0.0              0.0         0.00         0.00  180.00 48   4.48
 23 17  7  24.71  49.41   0.1844 29.3968   174.70  23 17  7       0.0              0.0         0.00         0.00    0.00 48   4.48
 23 13 13  24.71  49.41   0.1844 29.3968   174.70  23 13 13       0.0              0.0         0.00         0.00    0.00 24   4.48
 17 17 17  24.71  49.41   0.1844 29.3968   174.70  17 17 17       0.0              0.0         0.00         0.00  180.00  8   4.48
 20 18 12  24.72  49.44   0.1843 29.4307   174.90  20 18 12       0.0              0.0         0.00         0.00    0.00 48   4.48
 24 16  6  24.72  49.44   0.1843 29.4307   174.90  24 16  6       0.0              0.0         0.00         0.00    0.00 48   4.48
 22 18  8  24.78  49.57   0.1839 29.5663   175.71  22 18  8       0.0              0.0         0.01         0.00    0.00 48   4.45
 26 14  0  24.78  49.57   0.1839 29.5663   175.71  26 14  0       0.0              0.0         0.01         0.00    0.00 24   4.45
 24 14 10  24.78  49.57   0.1839 29.5663   175.71  24 14 10       0.0              0.0         0.01         0.00    0.00 48   4.45
 29  5  3  24.83  49.66   0.1836 29.6680   176.31  29  5  3       0.0              0.0         0.00         0.00  180.00 48   4.43
 23 15 11  24.83  49.66   0.1836 29.6680   176.31  23 15 11       0.0              0.0         0.00         0.00    0.00 48   4.43
 19 17 15  24.83  49.66   0.1836 29.6680   176.31  19 17 15       0.0              0.0         0.00         0.00  180.00 48   4.43
 27 11  5  24.83  49.66   0.1836 29.6680   176.31  27 11  5       0.0              0.0         0.00         0.00  180.00 48   4.43
 25 13  9  24.83  49.66   0.1836 29.6680   176.31  25 13  9       0.0              0.0         0.00         0.00    0.00 48   4.43
 25 15  5  24.83  49.66   0.1836 29.6680   176.31  25 15  5       0.0              0.0         0.00         0.00  180.00 48   4.43
 26 10 10  24.84  49.69   0.1835 29.7019   176.51  26 10 10       0.0              0.0         0.00         0.00    0.00 24   4.43
 26 14  2  24.84  49.69   0.1835 29.7019   176.51  26 14  2       0.0              0.0         0.00         0.00    0.00 48   4.43
 22 14 14  24.84  49.69   0.1835 29.7019   176.51  22 14 14       0.0              0.0         0.00         0.00    0.00 24   4.43
 21 19  9  24.95  49.90   0.1828 29.9393   177.93  21 19  9       0.0              0.0         0.00         0.00    0.00 48   4.39
 21 21  1  24.95  49.90   0.1828 29.9393   177.93  21 21  1       0.0              0.0         0.00         0.00  180.00 24   4.39
 26 12  8  24.96  49.93   0.1827 29.9732   178.13  26 12  8       0.0              0.0         0.00         0.00    0.00 48   4.38
 22 16 12  24.96  49.93   0.1827 29.9732   178.13  22 16 12       0.0              0.0         0.00         0.00    0.00 48   4.38
 28  8  6  24.96  49.93   0.1827 29.9732   178.13  28  8  6       0.0              0.0         0.00         0.00    0.00 48   4.38

@}

\end{document}
