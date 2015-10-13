
% nuweb formatted latex document for ID11 processing.
% 

\documentclass[11pt,notitlepage]{article}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage[a4paper]{geometry}

\begin{document}

\title{GraiTRE : GRAIn Translation REfinement}
\author{Jon Wright}
\date{Spring 2013}
\maketitle

\begin{abstract}
Try to fit the new diffractometer geometry with grains
that are not on the centre of rotation.
\end{abstract}

\tableofcontents

\section{Introduction}

After many attempts it seems to be quite difficult to get the 
refinement of diffractometer geometry correct using all the derivatives.
The purpose of this document is to figure things out step by
step and test each little piece separately prior to putting
it all together.
We would like the current approach to be compatible with the
previous work in ImageD11 and fable and also the current
instrumentation at the beamline.

\section{From the model to the data}

An approximate orientation matrix, unit cell and grain centre of
mass position are known. 
We want to refine the geometry.
Therefore we want to compute the positions of the diffraction 
spots on the detector, and their derivatives, and optimise
the various parameters.

The orientation and unit cell of the grain are described by
the conventional $(UB)$ matrix. 
This matrix is defined via the equation relating the 
$hkl$ indices of the diffraction peak to the reciprocal 
space scattering vector:

\begin{math}
 \begin{pmatrix} 
    g_0 \\ g_1 \\ g_2
 \end{pmatrix} 
 = 
   \begin{pmatrix} 
    ub_0 & ub_1 & ub_2 \\ 
    ub_3 & ub_4 & ub_5 \\ 
    ub_6 & ub_7 & ub_8 
    \end{pmatrix}  
   \begin{pmatrix} 
    h_0 \\ h_1 \\ h_2 
    \end{pmatrix}  
\end{math}
 
The subscripts correspond to the memory address locations as 
the vector is written in memory in the C programming language.
In fortran things will be the wrong way around, but we are
not using fortran, and everything here is in row major format.
That is to say matrices are stored by rows.
If we input the vector $(hkl)=(1,0,0)$ into that definition
we see that $ \begin{pmatrix} g_0 & g_1 & g_2 \end{pmatrix} =
\begin{pmatrix} ub_0 & ub_1 & ub_2 \end{pmatrix} $.
Thus a simple indexing algorithm involves picking three non-coplanar
difference vectors out of your data and writing down the 
orientation matrix $UB$ as rows. 


Habitually in the ImageD11 software we provided the UBI matrix, which
is the matrix inverse of UB. 
The reason is that ImageD11 was based on going mathematically from 
the data to the grains.

\begin{math}
 \begin{pmatrix} 
    h_0 \\ h_1 \\ h_2
 \end{pmatrix} 
 = 
   \begin{pmatrix} 
    ubi_0 & ubi_1 & ubi_2 \\ 
    ubi_3 & ubi_4 & ubi_5 \\ 
    ubi_6 & ubi_7 & ubi_8 
    \end{pmatrix}  
   \begin{pmatrix} 
    g_0 \\ g_1 \\ g_2 
    \end{pmatrix}  
\end{math}

You will notice there has been an interchange of row and column
vectors at this point.
In the ImageD11 code for indexing unknown cells this interchange
caused a great deal of confusion.
The concept of a row and column vectors as being ``different'' things
is difficult.

This matrix corresponds to the reciprocal space lattice vectors
stored as \emph{column} vectors. 
From http://en.wikipedia.org/wiki/Reciprocal\_lattice we paraphrase
to make $\mathbf{g_i}$ be the reciprocal lattice vectors and 
$\mathbf{a_i}$ be the real space unit cell vectors.

\[ \left[ \mathbf{g_0 g_1 g_2} \right]^T =   
   \left[ \mathbf{a_0 a_1 a_2} \right]^{-1} \]   

or, these are indeed rows as above, after the transpose:

\[ \begin{pmatrix} g_0 \\ g_1 \\ g_2 \end{pmatrix}=  
     \left[ \mathbf{a_0 a_1 a_2} \right]^{-1} \]   

or (TODO: this is confusing, add an example):

\[ \mathbf{ G }^T =   \mathbf{ A }^{-1} \]
\[ \left(\mathbf{ G }^T\right)^{-1} =   \left(\mathbf{ A }^{-1}\right)^{-1} \]
\[ \left( \mathbf{ G G^{-1} }\right)\mathbf{ G }^{-T} =   \mathbf{ A } \]
\[ \mathbf{ A } = \left( \mathbf{ G }\right) \left( \mathbf{G^{-1} G^{-T}}   \right)   \]
\[ \mathbf{ A } = \mathbf{ G } \left( \mathbf{G^{T} G}\right)^{-1}   \]

The row and column aspects remain confusing, as in software we currently
do not make a distinction between row and column vectors. 

We need to relate A and G to UB and UBI properly: there is an extra transpose.


\subsection{Finding the angle of diffraction}

Computation of the scan angle from the UB matrix.
This was quite a difficult problem in ImageD11. 
With hindsight it was a mistake not to have read the literature in 
more detail.
Seems the first nice derivation was due to Milch and Minor (1974). 
This is paraphrased later by D.J. Thomas, who covers the more
general case at the cost of introducing some unfamiliar notation.

We reproduce the story from Milch and Minor. The condition for diffraction is:

\[ \mathbf{| k_0 + R | = |k_0| } \]

That equation says the length of the incident scattering vector is
equal to the scattered vector, e.g., diffraction is elastic.
Here $\mathbf{k_0}$ is the incident wave-vector and $\mathbf{R}$ is the 
crystal lattice vector. Squaring both sides:
\[ \mathbf{| k_0 + R |^2 = |k_0|^2 } \]
\[ \mathbf{ k_0^2 + 2k_0R + R^2 - k_0^2} = 0  \]
\[ \mathbf{  R^2 + 2k_0R} = 0 \]
At this point Milch and Minor make some choice about the axis direction.
We turn to Thomas (1990, Goniometry) who has distinguished the two
different equations:
\[ \mathbf{  R^2 + 2k_0R} =  \mathbf{ R^2 + 2Rk_0}  \]
That is to say, the order of the dot product is reversed. 
The vector $\mathbf{R}$ is given by rotating the crystal lattice
vector according to the diffractometer rotation.
If it written as a row or coluwn vector, the rotation matrix is inversed, 
to act from the left on a column or from the right on a row vector. 
He distinguishes:

\begin{math}
 \begin{pmatrix} 
    R_0 \\ R_1 \\ R_2
 \end{pmatrix} 
 = 
    \Psi
   \begin{pmatrix} 
    g_0 \\ g_1 \\ g_2 
    \end{pmatrix}
\end{math} and \begin{math}
 \begin{pmatrix} 
    R_0 &  R_1 & R_2
 \end{pmatrix} 
 = 
   \begin{pmatrix} 
    g_0 & g_1 & g_2 
    \end{pmatrix}  
    \Psi^{-1}
\end{math}

Thus we get:

\[ \mathbf{  R^2 + 2k_0R} =  \mathbf{ R^2 + 2Rk_0} = 0  \]

\[ \begin{pmatrix} 
    R_0  &  R_1  &  R_2
 \end{pmatrix} 
 \begin{pmatrix} 
    R_0 \\ R_1 \\ R_2
 \end{pmatrix} 
+
2\begin{pmatrix} 
    R_0  &  R_1 &  R_2
 \end{pmatrix}
 \begin{pmatrix} 
    k_0 \\ k_1 \\ k_2
 \end{pmatrix} 
=
 \begin{pmatrix} 
    R_0  &  R_1  &  R_2
 \end{pmatrix} 
 \begin{pmatrix} 
    R_0 \\ R_1 \\ R_2
 \end{pmatrix} 
+
2\begin{pmatrix} 
    k_0  &  k_1 &  k_2
 \end{pmatrix}
 \begin{pmatrix} 
    R_0 \\ R_1 \\ R_2
 \end{pmatrix} 
\]

\[ \begin{pmatrix} 
    g_0 &  g_1 &  g_2 
    \end{pmatrix}  
    \Psi^{-1}
 \Psi
   \begin{pmatrix} 
    g_0 \\ g_1 \\ g_2 
    \end{pmatrix}
+
2\begin{pmatrix} 
    g_0  & g_1  & g_2 
    \end{pmatrix}  
    \Psi^{-1}
\begin{pmatrix} 
    k_0 \\ k_1 \\ k_2
 \end{pmatrix} 
= 0
\]
or,
\[ \begin{pmatrix} 
    g_0 &  g_1 &  g_2 
    \end{pmatrix}  
    \Psi^{-1}
 \Psi
   \begin{pmatrix} 
    g_0 \\ g_1 \\ g_2 
    \end{pmatrix}
+
2\begin{pmatrix} 
    k_0  &  k_1 &  k_2
 \end{pmatrix}
    \Psi
   \begin{pmatrix} 
    g_0 \\ g_1 \\ g_2 
    \end{pmatrix} =0
\]

\[ | \mathbf{ g} |^2 + 2 \begin{pmatrix} 
    g_0 &  g_1 &  g_2 
    \end{pmatrix}  
    \Psi^{-1}
\begin{pmatrix} 
    k_0 \\ k_1 \\ k_2
 \end{pmatrix} 
= | \mathbf{g} |^2
+
2\begin{pmatrix} 
    k_0  &  k_1 &  k_2
 \end{pmatrix}
    \Psi
   \begin{pmatrix} 
    g_0 \\ g_1 \\ g_2 
    \end{pmatrix} = 0
\]

This is a lot of symbols for apparently little gain so far.
One of the key points from Thomas' paper (1990, Goniometry)
is how to handle the rotation. 
He breaks it down into component parts and unifies the "axis-angle" 
view of rotations with the matrix form. 
Where we have written the operator $\Psi$ we normally mean it
as a 3x3 matrix. 
The insight is to decompose this operator into a sum of three different 3x3
matrices which correspond to the different projections of a vector
onto a rotation axis.
These three projections are the component along the axis, the component
in the plane of the vector and axis and the component making a 
right handed set. 
Then the rotation can be re-written as:

\[ \Psi = \Psi_{parallel} + \Psi_{plane}\cos\psi + \Psi_{perp}\sin\psi \]

\begin{figure}
\includegraphics[width=8cm]{decompose.png}
\caption{Reproduced from Thomas (1990, Goniometry)}
\end{figure}
For a rotation axis which has normalised direction $(l,m,n)$ these 
three matrices are:
\begin{math}
 \Psi_{p} = 
    \begin{pmatrix} 
        l^2  & ml  & nl \\
        lm  & m^2 &  nm \\
        ln  & mn &   n^2
    \end{pmatrix},
 \Psi_{c} = 
    \begin{pmatrix}
    1-l^2  & -ml &  -nl \\
    -lm  &  1-m^2 &  -nm \\
    -ln  &  -mn  &  1-n^2 
    \end{pmatrix},
 \Psi_{s} =
    \begin{pmatrix}
    0  & -n &  m \\
    n  & 0  & -l \\
    -m  & l &  0
    \end{pmatrix}
\end{math}

These are labelled with subscript $p$ for ``parallel'', $c$ for ``cosine'' component
and $s$ for ``sin'' component. Thus,

\begin{math}
 \Psi =
    \begin{pmatrix} 
        l^2  & ml  & nl \\
        lm  & m^2 &  nm \\
        ln  & mn &   n^2
    \end{pmatrix} +
    \begin{pmatrix}
    1-l^2  & -ml &  -nl \\
    -lm  &  1-m^2 &  -nm \\
    -ln  &  -mn  &  1-n^2 
    \end{pmatrix}
    \cos{\psi} + 
    \begin{pmatrix}
    0  & -n &  m \\
    n  & 0  & -l \\
    -m  & l &  0
    \end{pmatrix}
    \sin{\psi}
\end{math}

If we return to our problem of solving that equation we now find 
that everything is in terms matrices where we know all of the element
values and the only unknown is $\psi$. 
FIXME: what is k?
This solves the question of a general rotation axis. Amazing.
The derivative with respect to angle is also easier to find.

If we now consider the diffractometer as being a stack of rotations
we run into confusion.
The idea in ImageD11 was to consider only a single \emph{moving} axis.
Those rotations between the moving axis and the crystal should be
applied to the $\mathbf{UB}$ matrix prior to angle calculations.
Those between the moving axis and the floor are applied, in reverse,
to the incident beam direction. 
In this way we can always transform ourselves into a system where 
we can solve that equation, unless two axes are moving at the same time.
For the moment, we won't worry about that.

The solution is therefore found from:

\[ |\mathbf{g}|^2 + 2 \begin{pmatrix} k_0 & k_1 & k_2 \end{pmatrix}
\left[ \Psi_p + \Psi_c\cos{\psi} + \Psi_s \sin{\psi}
\right]
\begin{pmatrix} g_0 \\ g_1 \\ g_2 \end{pmatrix}
 = 0
\]

This is rougly the same thing as from gv\_general.py in ImageD11.
We can convert it to the form:
\[ A \sin{x} + B \cos{x} = C \sin{(x+\delta)} \]
where $ C = \sqrt{ A^2 + B^2 } $ and $\sin{\delta} = B/C $ and 
$\cos{\delta} = A/C$. 
This is most easily proved by substitution and expansion of the formulas,
google knows various pages that work through it.

Thomas goes further in (1992, Diffraction Geometry) in deriving a 
form of the final wavevector ($R$) in component form without resorting
to trig functions, but only a square root.



\subsection{Finding the point where a ray intersects a plane}
\label{rayplane}
Thomas (1992, Diffraction Geometry) invites us to consider the
matrix:

\[ \mathbf{d} = \left[ d^Y d^Z d^O \right] =
    \begin{pmatrix} 
            d_x^Y & d_x^Z & d_x^O \\
            d_y^Y & d_y^Z & d_y^O \\
            d_z^Y & d_z^Z & d_z^O
    \end{pmatrix} \]

The $Y$ and $Z$ correspond to the pixel directions in the detector
and the vector and $O$ is a vector from the origin (perhaps in the sample) 
to the point on the detector that we consider to be the origin.
Note that we do not require these vectors to be perpendicular, also
we do not need that $O$ is parallel to the beam.

The matrix of co-vectors is given by $\mathbf{dD = I}$:

\[ \mathbf{D} = \begin{pmatrix} ^YD \\ ^ZD \\ ^OD \end{pmatrix} =
    \begin{pmatrix} 
            ^YD_x & ^YD_y & ^YD_z \\
            ^ZD_x & ^ZD_y & ^ZD_z \\
            ^PD_x & ^PD_y & ^PD_z
    \end{pmatrix} \]

A scattered ray is given by:

\[ \begin{pmatrix} d^Y & d^Z & d^O \end{pmatrix}
   \begin{pmatrix} Q^Y \\ Q^Z \\ 1 \end{pmatrix} = \alpha \mathbf{R} \]

Here $\alpha$ is a constant of proportionality with dimensions of area. 
It is something like the conversion from millimeters to reciprocal
space angstroms.
Solving for $Q^Y$ and $Q^Z$ is then a simple matter of using $D$
and removing $\alpha$.


\[ \begin{pmatrix} Q^Y \\ Q^Z \\ 1 \end{pmatrix} = \alpha \mathbf{D.R} \]

\[ \begin{pmatrix} Q^Y \\ Q^Z \\ 1 \end{pmatrix} = \alpha 
    \begin{pmatrix} 
            ^YD_x & ^YD_y & ^YD_z \\
            ^ZD_x & ^ZD_y & ^ZD_z \\
            ^PD_x & ^PD_y & ^PD_z
    \end{pmatrix}
    \begin{pmatrix} R_x \\ R_y \\ R_z  \end{pmatrix}
 \]

This expands to give 3 equations for each of our rows:

\[ Q^Y = \alpha \left( {^Y}D_x R_x + {^Y}D_y R_y + {^Y}D_z R_z \right) \]
\[ Q^Z = \alpha \left( {^Z}D_x R_x + {^Z}D_y R_y + {^Z}D_z R_z \right)\]
\[ Q^P = \alpha \left( {^P}D_x R_x + {^P}D_y R_y + {^P}D_z R_z \right) = 1 \]

We use the last one to solve for alpha and back substitute that into
the others:

\[ Q^I = \frac{ {^I}D_x R_x + {^I}D_y R_y + {^I}D_z R_z } 
              { {^P}D_x R_x + {^P}D_y R_y + {^P}D_z R_z }  \]
...for $I$ is $Y$ or $Z$.

Cool.


\section{Computing from the data towards the model}

This was the approach taken in ImageD11 in the past. 
It seems to be relatively straightforward to keep that code
compatible with anything new here.

\section{Combining translations with the other parameters}

The laboratory co-ordinates can be computed from the detector position.

The grain position is computed from the translation.

We can consider the vector from grain to detector in either the frame
of reference of the detector or the diffractometer.

In the detector frame we have the vector:

\[ dX = L + \Omega T \]

where $dX$ is the vector from grain to detector, $L$ are the lab co-ordinates
on the detector plane and $\Omega$ is the rotation of the diffractometer.

The scattering vector, $k$, is 
\[ k = \left( dX/|dX| - s_0 \right) / \lambda \]

The g-vector, in the crystal frame, is then rotated back:

\[ g = \Omega^{-1} k = \Omega^{-1} \left( dX/|dX| - s_0 \right) / \lambda \]


\[ \lambda g = \Omega^{-1}dX/|dX| - \Omega^{-1} s_0  \]

We substitute for $dX$ and $g = UB.h$:

\[ \lambda g = \lambda UB.h = \Omega^{-1}( L + \Omega T ) /|dX| - \Omega^{-1} s_0  \]

\[ \lambda |dX| UB.h = \Omega^{-1} L + T - \Omega^{-1} |dX| s_0  \]

\[ T = \lambda |dX| UB.h - \Omega^{-1} L + \Omega^{-1} |dX| s_0  \]

So we see there is a scale factor, per reflection, relating laboratory
co-ordinates to reciprocal angstroms, $ \lambda |dX| $.

This equations suggests we can compute $T$ in a single step as the reciprocal
space origin, possibly cycling a couple of times to take account of the effect
of $T$ on $dX$.

\section{Putting this into code}

Our aim is to fit our diffractometer geometry.
This can be in terms of the extracted peak positions from a 
peak search of a stack of images.
This can also be in terms of a point spread function and the 
individual pixel intensities in the images.

Eventually we need our code to be robust and well tested and 
reasonably efficient. 
We would like to be able to call it on a graphics card inside
an OpenCL kernel, from the C language and from python.

\subsection{Experiment descriptions}

In fable the geometry is pretty much hard wired into the software.

In programs like XDS a more general description is available.

We want to write a fable (ImageD11) geometry description in a way 
that will work with the old code and allow new instrument definitions
and also be pretty fast.

\subsubsection{Detectors}

The output of peaksearch is the peak positions in the pixel slow and fast
directions and also the observed omega angles. 
The instrument ``parameters'' are a flip matrix, the beam center (2) in pixel
co-ordinates, the pixel sizes (2), the detector tilts (3) and the 
sample to detector distance (1).
This makes a total of 8 parameters with the problem of things blowing up
if the detector is on a two theta arm.

Some equivalent description could be made, perhaps adding in a detector
two theta arm (horizontal or vertical).

Eventually we should write these parameters into a matrix description 
as in section \ref{rayplane}.

For refinement it seems we need something like a parameter list,
a function converting that list to a single matrix and
a another function giving derivatives of the matrix elements
with respect to each parameter.







Rotation axis requires minimum of 3 direction cosines.

Build a stack of rotations to represent diffractometer.

Apply rotations to incident beam and UB to reduce to single
axis form.

\section{ Coding in python }

Try to write an initial version in python because
we are most fluent in this language at the moment
and it should be the easiest to debug (interactive 
prompt).

So, what is the overall computation? 
We want to compute the (sc, fc, angle) for a series
of indexed diffraction spots in terms of a model.
That is detector position and scan rotation angle.

We will start with just the angle part. 
The model of the diffractometer is a series of rotation
axes stacked on top of each other. 
Represent this as a list of axis directions.
One axis is noted as the moving axis.
We also input the wavelength, UB matrix and hkl indices:

@o compute_omega.py @{

"""
Module to implement the maths for the grain translation refinement
"""

import numpy as np
import math

def rodrigues_formula( vector, axis, angle):
    """
    vector = [3,] thing to be rotated
    axis   = normalised vector along axis
    angle  = angle to rotate by (degrees)

    Implements the Rodrigues Rotation Formula from
    wikipedia: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    v_rot = v cos(theta) + (k x v) sin(theta) + k(k.v)(1 - cos(theta))
    """
    cost = math.cos(math.radians(angle))
    sint = math.sin(math.radians(angle))
    nhat = np.asarray(axis) / math.sqrt( np.dot(axis, axis) )
    kxv = np.cross(nhat, vector)
    kkdv = nhat * np.dot(nhat, vector)
    vrot = cost*np.asarray(vector) + kxv*sint + kkdv*(1-cost)
    return vrot
    
def rotation_matrix( axis, angle ):
    """
    Convert axis angle to rotation matrix
    """
    return np.array( [rodrigues_formula(v, axis, angle) for v in np.eye(3)] )

def get_angle( scanaxis, axes, axpos, incident_beam, gvec ):
    """
    axis = integer pointer into list of axis
    axes = vectors along axis directions
    a1   = positions of non-scanning axes
    s0   = incident beam wave vector
    g0   = reciprocal lattice vector at a0
    """
    assert scanaxis > 0 and scanaxis < len(axes) # in range
    matrices = [ rotation_matrix( axis, angle ) for (axis, angle) in
    	     zip( axes, axpos ) ]
    # Rotate the incident beam according to the first axes
    # as the list starts at the floor i increments
    rot_inbeam = incident_beam
    for i in range(0, scanaxis):
        rot_inbeam = np.dot(matrices[i], rot_inbeam)
    # Rotate crystal lattice vector to scan axis co-ords
    # as the list ends at the crystal i decrements
    rot_gvec = gvec
    for i in range(len(axes), scanaxis ):
        assert i != scanaxis
        rot_gvec = np.dot( matrices[i], rot_gvec)
    # We are ready to solve a single axis angle problem now
    angle = axis_solve( axes[scanaxis], rot_inbeam, rot_gvec )


def axis_solve( axis, s0vec, g0vec ):
    """
    axis = rotation axis
    s0vec = incident beam direction
    g0vec = crystal lattice at angle = 0 
    """
    return np.dot( np.dot(axis, s0vec), g0vec)
    
if __name__ == "__main__":
    print  rodrigues_formula( (0, 0, 1), 
    	   		     (0, 1, 0), 
                             30)
    print  rotation_matrix( (0, 1, 0), 30)


@}


\subsubsection{ Testing the code }

We try to set up a series of test cases to verify whether
out computations are correct and meaningful.


@o test_compute_omega.py 
@{

"""
Test cases for the compute_omega module
"""

import unittest, numpy.testing as npt
import compute_omega as co

# DELTA is the smallest value which can be added to
# one without upsetting it. e.g. the precision of a
# double precision number
DELTA = 1.
while (1. + DELTA) != 1.: 
    DELTA /= 2


class TestComputeOmega(unittest.TestCase):
    """
    Test case class
    """      
    def test_rodrigues_formula(self):
    	""" 
	Testing rotation of vectors about vectors:
      	rodrigues_formula(vector, axis, angle)
	"""
	xvec = (1, 0, 0) 
	yvec = (0, 1, 0) 
	zvec = (0, 0, 1) 
	# X -> Y on +90 around Z
	vcalc = co.rodrigues_formula(xvec, zvec, 90)
	npt.assert_allclose( vcalc , yvec , atol = DELTA*2 )
	# Z -> X on +90 around Y
	vcalc = co.rodrigues_formula(zvec, yvec, 90)
	npt.assert_allclose( vcalc , xvec , atol = DELTA*2 )
	# Y -> Z at +90 around X
	vcalc = co.rodrigues_formula(yvec, xvec, 90)
	npt.assert_allclose( vcalc , zvec , atol = DELTA*2 )
	# Y -> -Y at +180 around X
	vcalc = co.rodrigues_formula(yvec, xvec, 180)
	npt.assert_allclose( -vcalc , yvec , atol = DELTA*2 )
	# Y -> -Z at +270 around X
	vcalc = co.rodrigues_formula(yvec, xvec, 270)
	npt.assert_allclose( -vcalc , zvec , atol = DELTA*2 )
	



if __name__ == '__main__':
    unittest.main()	  

@}


\section{Examples}

\section{Conclusion}

\section{Building the code}

@o bld.bat
@{
nuweb graitre
pdflatex graitre
bibtex graitre
pdflatex graitre
pdflatex graitre
@}

@o makefile -t
@{
 
all: graitre.pdf pylint.out test.out

graitre.pdf: graitre.w
	nuweb graitre
	pdflatex graitre


pylint.out: graitre.w
	nuweb graitre
	pylint compute_omega.py | tee pylint.out
	pylint test_compute_omega.py | tee -a pylint.out

test.out: graitre.w
	nuweb graitre
	python test_compute_omega.py | tee tests.out

@}



\nocite{*}
\bibliographystyle{plain}
\bibliography{graitre}


\end{document}
