% nuweb formatted latex document for ID11 processing.
% 

\documentclass[11pt,notitlepage]{article}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage[a4paper]{geometry}

\begin{document}

\title{WWM : Wafer Wavelength Monitor}
\author{Jon Wright}
\date{Summer 2013}
\maketitle

\begin{abstract}
Data processing code for the wafer wavelength monitor scan.
\end{abstract}

\tableofcontents

\section{Introduction}

A fast scan has been developed to measure the extinction diffraction
of a silicon wafer in order to make a wavelength monitor for ID11.
Using the ESRF MUSST card data are recorded for the incident and 
transmitted beam intensity and these are written to an ascii spec file.

In order to extract the wavelength, crystal orientation and axis direction
from the data some processing is required. 

The following steps are identified:

Normalisation of the scan data

Extraction of peak positions, widths and heights as a function of angle

Assignment of hkl indices to the peaks

Refinement of the geometry (3 orientation, one axis tilt, one wavelength)

\section{Theory}

There is a nice derivation in Milch and Minor (1974).
This is also something like the Bond method (ref).

We reproduce the Milch and Minor story
% FIXME - get the notation to match fable/ImageD11.
. The condition for diffraction is:

\[ \mathbf{| k_0 + R | = |k_0| } \]

Here $\mathbf{k_0}$ is the incident wave-vector and $\mathbf{R}$ is a 
crystal lattice vector. 

\[ \mathbf{| k_0 + R |^2 = |k_0|^2 } \]
\[ \mathbf{ k_0^2 + 2k_0R + R^2 - k_0^2} = 0  \]
\[ \mathbf{  R^2 + 2k_0R} = 0 \]

Once the hkl is assigned to the peak then $|\mathbf{R}|^2$ is a constant
which does not depend on rotations. 
During the experiment the component of the vector along the rotation axis,
$k_z$, does not change, so that the rotated vector $R_z$ is also constant.
We define a co-ordinate system where the rotation axis defines z and the 
beam is incoming in the x-z plane.
Ideally, if the beam is perpendicular to the axis, the it will be along x.
Otherwise we get $k_0 = k_x  + k_z= x.\cos{B}+y.\sin{B}$, where B is the angle
between the beam and axis, hopefully rather small.
If we pop that into the eqn above we get:

\[ |R|^2 + (k_x + k_z).(R_x + R_y + R_z) = 0 \]

\[ |R|^2 + k_x.R_x + k_z.R_z = 0 \]

\[ k_x.R_x = - |R|^2 - k_z.R_z \]

\[ R_x = - \frac{|R|^2 - k_z.R_z }{k_x} \]

If the orientation is known, at this point we know all those terms.

The component along the axis is constant, $R_z$, and so we can get the
$R_y$ component from $R_y=+/0 \sqrt{|R|^2}-R_x^2-R_z^2$.

We have two positions of $\mathbf{R}$ where the diffraction can happen.
We solve for the rotation angle, $\phi$, since:

\[ R_y = R_{0y} \cos{\phi} - R_{0x} \sin{\phi} \]

\subsection{Extinction diffraction}

We are concerned by the overall length of the diffraction vector and
the projection of it onto the rotation axis.
For a single peak we can observe both $g$ and $-g$.
From these four observations we will find peaks at $\phi$ and
$\phi+180$ if the axis is perpendicular to the beam.
If this is not the case there is an offset.
The difference in angle between the two Friedel pairs gives us a single 
observation to constrain the length and axis projection (two unknowns).

If we tweak the wavelength we change $|R|^2$ but keep the same direction
cosines. 
This tells us the sign of the solution in the square root equation as the
peaks shift either to lower or higher angles. 

Example data from 65 and 65.1keV show this in action. 
Peaks shift by their width (about 100 eV) either to positive or negative.
Looks like some Freidel pairs get closer together in angle, some further apart. 
None of them look like staying where they are...

\section{Axis alignment procedure}

If a tilt is installed under the axis measure a pair of sharp peaks at 0 and 180 
degrees.
Plot the separation of these peaks (-180) versus the tilt angle.
When the axis is perpendicular to the beam the peaks will be 180 degrees apart.








\section{spec macros}

Include the spec macro file

\section{MUSST program}

Include the musst program

\section{Building the code}

@o bld.bat
@{
nuweb wwm
pdflatex wwm
@}

\section{ Literature }
The Indexing of Single-Crystal X-ray Rotation Photographs.
J. R. Milch and T. C. Minor.
J. Appl. Cryst. (1974) 7, 502

\end{document}
