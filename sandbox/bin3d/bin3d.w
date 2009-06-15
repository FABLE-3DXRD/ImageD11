

\documentclass[10pt,a4paper,notitlepage]{article}
\usepackage{graphics} % For the pictures
\usepackage{anysize}  % Try to circumvent Latex default margins

\newcommand{\var}[1]{\textbf{\textsf{#1}}} % highlight variables in text
\newcommand{\code}[1]{\textbf{\textsf{#1}}} % highlight code in text
\newcommand{\param}[1]{\textbf{\textsf{#1}}} % ... parameters ...
\begin{document}

\marginsize{2cm}{2cm}{2cm}{2cm} % Needs anysize
\pagestyle{headings}            % These are ugly - fix them somehow?
\pagenumbering{arabic}          % Page numbers

\title{Rebinning in 2D and 3D}
\author{JPW}
\date{Summer 2009}
\maketitle

\begin{abstract}
Some code for cutting up pixels in 2D and 3D.
\end{abstract}

\tableofcontents

% \newpage    % Only need this if contents don't run over first page

\section{Introduction}
We have a rebinning program for 1D which is used on ID31. 
It splits up a 1D pixel onto a grid of regular bins.

For doing radial integration and also projection diffraction spot shapes 
into sample co-ordinates it would be useful to have a 2D rebinning routine. 
This should take a 2D polygon and cut it up into little pieces which 
exactly match the 2D rebinning target grid. 
The polygon will be represented as a list of connected edges.
Each edge is a 1D rebinning problem (it seems)

When projecting images into a sample volume or reciprocal space volume 
it could be interesting to make a 3D overlap computation. 
In this case the 3D polyhedron should be represented as a collection of faces.
The rebinning is then a question of cutting up into smaller polyhedra. 
It would seem that the planes representing the regular grid in 3D 
will cut the polyhedral faces to give a set of 2D rebinning problems.

What do we notice? It looks like an N dimensional object is representable 
as a list of N-1 dimensional objects which describe the boundary. 
Cutting up an N dimensional object can be done by manipulating the (N-1) 
dimensional things recursively until we reach lines and points which are 
the base case. However, 3D is probably enough.

\section{Technical requirements}
    
We'd like routines which are easy to call and re-use. 
Also well documented and easy to understand/remember.
The code should be able to run as fast as the computing hardware allows.
It should not fuck about with time.
The results should be demonstrably correct. 
Testcases and proofs are required.
Ideally, it will run on GPU platforms just as happily as CPU.

\section{Use cases}

Input is a 2D data image

Output is a 2D array which will receive the rebinned data [bounds big enough or not flag?]
The mapping from one array to the other.
* Function to call with pixel co-ordinates
* Pixel by pixel tabulated
* 2D affine matrix for regular splats

Input is a series of 2D images with "box" pointspread

Output is a high resolution reconstruction of 1D radial integrated from doing ART on the 1D integration

Input is a 3D image series

Output is a 3D reciprocal space volume

Input is a single 2D image

Output is the pixel by pixel projections into a 3D volume to supply to an ART program.

Input is a 1D list of counts

Output is a 1D list of counts
Mapping is change from d-spacing to two theta


\section{Plan}

Inputs:
source dimensions and geometry
Mapping function (dimensions, i, j, k ,...) -> ??? How to define the function and grid ???
destination dimensions and geometry

Outputs:
For each input pixel -> destination and multiplier fractions

Execution
Optionally a data array can be rebinned on-the-fly computing geometry as we go. This avoids saving potentially large intermediate data.
Alternatively we can save the "plan" in a file and execute it later without computing the geometry again.

\section{Implementation}

A difficult question. Historically Jon has tried Fortran77, 90, and python/C. 
Learning from this experience, fortran is just too hard to port/call/compile.
Python is too hard to call/distribute.
This leaves the C language for the "interesting" parts of the code.
Actual driver programs can be in any language that calls C. 
So what are the functions we need?

\subsection{Mappings 1D}

The functions to define the rebinning co-ordinate transform. 
In the 1D case we had a set of functions:


@o oneD.h
@{

/** 
 * A simple linear 1D mapping
 */
typedef struct  {
    float lowLimit ;
    float highLimit;
    float stepSize;
    } oneD ;

 int ibin( oneD*, float);

/**
 * the upper boundary for this bin == binl(n+1) 
 */
 float binh( oneD*, int);

/**
 * the lower boundarry for this bin ==binh(n-1) 
 */
 float binl( oneD*, int);

/**
 * the bin centre == 0.5(binl+binh) 
 */
 float bincen( oneD*, int);

@}

These should be relatively straightforward to implement. 
We will need some mapping data defined somewhere.

@O oneD.c
@{
#include <math.h>
#include "oneD.h"
/**
 * the bin number for this point 
 */
int ibin( oneD* mapping, float x ){
 if (  (x < mapping->lowLimit ) || (x > mapping->highLimit ) ){
    return -1;
    }
 return floor( 0.5 + (x - mapping->lowLimit ) / mapping->stepSize );
}

@}


\section{Building and testing this document}

We'll start by making a makefile for the whole thing.

@O Makefile -t
@{
CC = gcc
CFLAGS = -c -Wall 
AR = ar rvu
RANLIB = ranlib

SOURCES = oneD.c 
OBJECTS =  $(SOURCES:.c=.o)

all: bin3d.dvi libbin3d.a

bin3d.dvi: bin3d.tex
	latex bin3d.tex

bin3d.tex: bin3d.w
	nuweb bin3d.w

libbin3d.a : $(OBJECTS)
	$(AR) $@@ $(OBJECTS)
	$(RANLIB) $@@

.c.o:
	$(CC) $(CFLAGS) $< -o $@@


@}

\subsection{Testing the C code}

... need a C test case infrastructure ...

\enddocument

