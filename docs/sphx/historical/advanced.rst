========
Advanced
========

Take advanced to mean without graphical interface, without prejudice to people
who prefer touching mice to keyboards. The program was mainly developed for
online analysis during experiments to get an idea if the setup and sample etc
are OK. Other programs exist which might to a better job, and you should
consider trying them for serious data analysis. Notably, `Grainspotter
<https://sourceforge.net/p/fable/wiki/grainspotter/>`_ reads ImageD11 gve files
and tends to make a better indexing. Nevertheless, the free nature of ImageD11
and convenience of the python language could tempt you to process large volumes
of data using the functions in the program. Currently the indexing module is the
most easily accessible for script writers. The other modules should be separated
from the gui very soon.

Using the indexing module interactively
---------------------------------------

An example interactive session is given here. A trial orientation is generated
from one g-vector file and then used to index another::

  fable4:~/id11/lab6grains % python
  Python 2.4 (#3, Mar 11 2005, 10:10:40)
  [GCC 3.3 20030226 (prerelease) (SuSE Linux)] on linux2
  Type "help", "copyright", "credits" or "license" for more information.
  >>> from ImageD11 import indexing
  >>> myindexingobject = indexing.indexer()
  >>> dir(myindexingobject)   # Gives a list of commands/properties
  ['__doc__', '__init__', '__module__', 'assigntorings', 'cosine_tol',
  'coverage', 'ds_tol', 'find', 'friedelpairs', 'getind', 'gv', 'hkl_tol',
  'max_grains', 'minpks', 'readgvfile', 'refine', 'ring_1', 'ring_2',
  'saveindexing', 'saveubis', 'score', 'scores', 'scorethem', 'ubis',
  'uniqueness', 'unitcell', 'wavelength']
  >>> # Read in some filtered peaks (only the strong ones)
  >>> myindexingobject.readgvfile("filt.gve")
  Got wavelength from gv file of  0.3333
  Read your gv file containing (112, 3)
  >>> myindexingobject.assigntorings()
  Maximum d-spacing considered 0.991615
  Ring assignment array shape (112,)
  Ring     (  h,  k,  l) Mult  total indexed to_index
  Ring 0   ( -1,  0,  0)    6     11      0     11
  Ring 1   ( -1, -1,  0)   12     39      0     39
  Ring 2   ( -1, -1, -1)    8     17      0     17
  Ring 3   ( -2,  0,  0)    6     13      0     13
  Ring 4   ( -2, -1,  0)   24     13      0     13
  Ring 5   ( -2, -1, -1)   24      4      0      4
  Ring 6   ( -2, -2,  0)   12      2      0      2
  Ring 7   ( -3,  0,  0)   30      8      0      8
  Ring 8   ( -3, -1,  0)   24      3      0      3
  Ring 9   ( -3, -1, -1)   24      1      0      1
  Ring 10  ( -2, -2, -2)    8      0      0      0
  Ring 11  ( -3, -2,  0)   24      0      0      0
  Ring 12  ( -3, -2, -1)   48      0      0      0
  Ring 13  ( -4,  0,  0)    6      0      0      0
  Using only those peaks which are assigned to rings for scoring trial matrices
  Shape of scoring matrix (111, 3)
  >>> # Choose which rings to use for generating trial orientations
  >>> print myindexingobject.ring_1
  1
  >>> print myindexingobject.ring_2
  2
  >>> myindexingobject.ring_1=0
  >>> myindexingobject.ring_2=0
  >>> print myindexingobject.ring_1, myindexingobject.ring_2
  0 0
  >>> # Find the trial orientations or grains
  >>> myindexingobject.find()
  hkls of rings being used for indexing
  Ring 1: [(-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
  Ring 2: [(-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
  Possible angles and cosines between peaks in rings:
  90.0 6.12303176911e-17
  Number of peaks in ring 1: 11
  Number of peaks in ring 2: 11
  Minimum number of peaks to identify a grain 10
  Percent done 90.909%   ... potential hits 10
  Number of trial orientations generated 10
  Time taken 0.00019097328186
  >>> # Test those orientations against the data
  >>> myindexingobject.scorethem()
  Scoring 10 potential orientations
  Tested        9    Found        2     Rejected        0 as not being unique
  Number of orientations with more than 10 peaks is 2
  Time taken 0.00936889648438
  UBI for best fitting
  [[ 2.58780454 -3.00470196  1.24755448]
   [-1.09995311  0.69776998  3.94961316]
   [-3.05800078 -2.78861403 -0.35888152]]
  Indexes 42 peaks, with <drlv2>= 0.00370687876058
  That was the best thing I found so far
  Number of peaks assigned to rings but not indexed =  11
  >>> # Look at the orientations found:
  >>> for u in myindexingobject.ubis: print u
  ...
  [[ 2.59661965 -2.99631727  1.24883939]
   [-1.09621729  0.69589138  3.94892308]
   [-3.05545818 -2.79602729 -0.35546775]]
  [[ 0.16666909 -3.31769292  2.49901696]
   [-1.36226123  2.31902622  3.16959182]
   [-3.92382675 -0.94603473 -0.9942598 ]]
  >>> # see what the refined cell parameters are
  >>> indexing.ubitocellpars(myindexingobject.ubis[0])
  (4.156916, 4.156915999, 4.156916, 90.0, 89.9999,89.9999)
  >>> indexing.ubitocellpars(myindexingobject.ubis[1])
  (4.156916, 4.15691600, 4.1569159,90.00, 89.9999, 89.9999)
  >>> # These look a bit too good to be true...?
  
Got an orientation now - test it on another dataset (could be another 
rotation or a different detector etc)::

  >>> # Anyway, apply these to another dataset with weak peaks in it too (peak.gve)...
  >>> another_indexing_object = indexing.indexer()
  >>> another_indexing_object.readgvfile("peak.gve")
  Got wavelength from gv file of  0.3333
  Read your gv file containing (1734, 3)
  >>> another_indexing_object.assigntorings()
  Maximum d-spacing considered 1.177329
  Ring assignment array shape (1734,)
  Ring     (  h,  k,  l) Mult  total indexed to_index
  Ring 0   ( -1,  0,  0)    6     71      0     71
  Ring 1   ( -1, -1,  0)   12    160      0    160
  Ring 2   ( -1, -1, -1)    8     88      0     88
  Ring 3   ( -2,  0,  0)    6     65      0     65
  Ring 4   ( -2, -1,  0)   24    203      0    203
  Ring 5   ( -2, -1, -1)   24    143      0    143
  Ring 6   ( -2, -2,  0)   12     76      0     76
  Ring 7   ( -3,  0,  0)   30    191      0    191
  Ring 8   ( -3, -1,  0)   24    149      0    149
  Ring 9   ( -3, -1, -1)   24    127      0    127
  Ring 10  ( -2, -2, -2)    8     25      0     25
  Ring 11  ( -3, -2,  0)   24     95      0     95
  Ring 12  ( -3, -2, -1)   48    142      0    142
  Ring 13  ( -4,  0,  0)    6     13      0     13
  Ring 14  ( -4, -1,  0)   48     67      0     67
  Ring 15  ( -4, -1, -1)   36     47      0     47
  Ring 16  ( -3, -3, -1)   24     17      0     17
  Ring 17  ( -4, -2,  0)   24     13      0     13
  Ring 18  ( -4, -2, -1)   48     27      0     27
  Ring 19  ( -3, -3, -2)   24      5      0      5
  Using only those peaks which are assigned to rings for scoring trial matrices
  Shape of scoring matrix (1724, 3)
  >>> another_indexing_object.hkl_tol=0.1 # Make a bigger tolerance in case grain rotates
  >>> test_orientation = myindexingobject.ubis[0]
  >>> print test_orientation
  [[ 2.59661965 -2.99631727  1.24883939]
   [-1.09621729  0.69589138  3.94892308]
   [-3.05545818 -2.79602729 -0.35546775]]
  >>> new_orientation = another_indexing_object.refine(test_orientation)
  >>> print new_orientation
  [[ 2.58919314 -3.00117096  1.24689206]
   [-1.09852789  0.69658936  3.95070245]
   [-3.05916685 -2.79050878 -0.3580728 ]]
  >>> indexing.ubitocellpars(new_orientation)
  (4.1552001238, 4.1593328915, 4.15615895,89.993, 89.975, 90.029)
  >>> # These look like they were really refined. Maybe some peaks are
  >>> # used that should not be due to the hkl_tol of 0.1 (a bit big)
  >>> indexed_peaks = another_indexing_object.getind(new_orientation,tol=0.05)
  >>> print indexed_peaks.shape
  (339,)
  >>> # So 339 peaks seem to be indexed with the peak.gve dataset and tol=0.05
  >>> # Have a look at the hkls generated by the UBI matrix
  >>> for k in indexed_peaks:
  ...    print matrixmultiply(new_orientation,another_indexing_object.gv[k])
  ...
  [ 0.10168799 -0.9761811   0.13347752]
  [  9.13527496e-04   1.97985539e-04  -9.97830257e-01]
  [  9.97934113e-01   7.92302450e-04   3.49838267e-03]
  [ 0.99452548 -0.06653118 -0.05044758]
  [ -3.98020413e-04  -9.99843377e-01   1.47886732e-04]
  ...<many lines snipped>..
  >>> # etc
  

Transformations and peak merging without the gui
================================================

Use the history from Help>History to capture the commands run when you click on
the tk gui.

