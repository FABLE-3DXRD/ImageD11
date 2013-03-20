.. ImageD11 documentation master file, created by
   sphinx-quickstart on Thu Sep 13 12:28:22 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ImageD11's documentation!
====================================

ImageD11 is for identifying individual grains in spotty area detector diffraction images. 
It finds the spots and helps you to try to index them. 
After you index the spots you can then work on extracting useful information on a grain 
by grain basis. 
ImageD11 is intended to be complementary to the rest of the FABLE software. ImageD11 has 
a strong focus on fast preliminary data analysis which can be carried 
out during experiments in order to make decisions about the experimental setup etc.

The package includes facilities for peak searching in two dimensional 
area detector images using a simple threshold. 
A graphical interface is available to assist with the calibration 
of the extracted peak positions and their transformation into 
scattering vectors (or g-vectors) in reciprocal space.
In cases where "relatively few" grains are illuminated by the beam the 
package is capable of indexing the diffraction spots and refining 
unit cell parameters and orientation matrices. 
ImageD11 is mainly written using the python language, with a few 
extension modules in c for increased performance. 



Contents:

.. toctree::
   :maxdepth: 1

   installation
   programs
   peaksearching
   postprocessing
   filtering
   calibration
   indexing
   index_unknown
   refinement
   advanced
   fileformats
   batchprocessing


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

