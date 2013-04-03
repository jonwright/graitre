
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
orientation matrix. 

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

You will notice there has been and interchange of row and column
vectors at this point.
In the ImageD11 code for indexing unknown cells this interchange
caused a great deal of confusion.
The concept of a row and column vectors as being "different" things
is difficult.

This matrix corresponds to the reciprocal space lattice vectors
stored as \emph{column} vectors. 
From http://en.wikipedia.org/wiki/Reciprocal\_lattice we paraphrase
to make $\mathbf{g_i}$ be the reciprocal lattice vectors and 
$\mathbf{a_i}$ be the real space unit cell vectors.

\[ \left[ \mathbf{g_0 g_1 g_2} \right]^T =   
   \left[ \mathbf{a_0 a_1 a_2} \right]^{-1} \]   

or:

\[ \mathbf{ G }^T =   \mathbf{ A }^{-1} \]
\[ \left(\mathbf{ G }^T\right)^{-1} =   \left(\mathbf{ A }^{-1}\right)^{-1} \]
\[ \left( \mathbf{ G G^{-1} }\right)\mathbf{ G }^{-T} =   \mathbf{ A } \]
\[ \mathbf{ A } = \left( \mathbf{ G }\right) \left( \mathbf{G^{-1} G^{-T}}   \right)   \]
\[ \mathbf{ A } = \mathbf{ G } \left( \mathbf{G^{T} G}\right)^{-1}   \]

The row and column aspects remain confusing, as in software we currently
do not make a distinction between row and column vectors. 

We need to relate A and G to UB and UBI properly.


\subsection{Finding the angle of diffraction}

Computation of the scan angle from the UB matrix.
This was quite a difficult problem in ImageD11. 
With hindsight it was a mistake not to have read the literature in 
more detail.
Seems the first nice derivation was due to Milch and Minor (1974). 
This is paraphrased later by D.J. Thomas, who covers the quite
general case at the cost of introducing some unfamiliar notation.

We reproduce the story. The condition for diffraction is:

\[ \mathbf{| k_0 + R | = |k_0| } \]

Here $\mathbf{k_0}$ is the incident wave-vector and $\mathbf{R}$ is the 
crystal lattice vector. 
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

\includegraphics[width=8cm]{decompose}

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

These are labelled with subscript $p$ for "parallel", $c$ for cosine component
and $s$ for sin component. Thus,

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
This solves the question of a general rotation axis. Amazing.
The derivative with respect to angle is also made abundantly clear.

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
\[ A \sin{x} + B \cos{x} = C \sin{x+\delta} \]
where $ C = \sqrt{ A^2 + B^2 } $ and $\sin{\delta} = B/C $ and 
$\cos{\delta} = A/C$. 
This is most easily proved by substitution and expansion of the formulas,
google knows various pages that work through it.

Thomas goes further in (1992, Diffraction Geometry) in deriving a 
form of the final wavevector in component form without resorting
to trig functions, but only a square root.




\subsection{Finding the point where a ray intersects a plane}

\section{Computing from the data towards the model}



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



\nocite{*}
\bibliographystyle{plain}
\bibliography{graitre}


\end{document}
