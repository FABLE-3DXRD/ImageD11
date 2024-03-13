.. ImageD11 documentation master file, created by
   sphinx-quickstart on Thu Sep 13 12:28:22 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================================
Welcome to ImageD11's documentation!
====================================

Beware: you might find some of this information is out of date or wrong. Please
fix anything you find that is wrong and send a pull request. This documentation
lives in the docs/sphx part of the sources.

ImageD11 is for identifying individual grains in spotty area detector
diffraction images. It finds the spots and helps you to try to index them. After
you index the spots you can then work on extracting useful information on a
grain by grain basis. ImageD11 is intended to be complementary to the rest of
the FABLE software. ImageD11 has a strong focus on fast preliminary data
analysis which can be carried out during experiments in order to make decisions
about the experimental setup etc.

The package includes facilities for peak searching in two dimensional area
detector images using a simple threshold. A graphical interface is available to
assist with the calibration of the extracted peak positions and their
transformation into scattering vectors (or g-vectors) in reciprocal space.In
cases where "relatively few" grains are illuminated by the beam the package is
capable of indexing the diffraction spots and refining unit cell parameters and
orientation matrices. ImageD11 is mainly written using the python language, with
a few extension modules in c for increased performance.


Here is a pdf file with a worked example for making a centre of mass grain map
:download:`com_guide.pdf<../com_guide.pdf>`


Contents:

.. toctree::
   :maxdepth: 1

   installation
   programs
   peaksearching
   filtering
   calibration
   indexing
   index_unknown
   refinement
   advanced
   parallel
   fileformats
   batchprocessing
   rsv_mapper
   related
   changelog


Indices and tables
==================

All the code docstrings can be found via:

.. toctree::
   :maxdepth: 1

   api/modules

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
