index_unknown.py
================
Warning: This code is developmental, so it might not work for you, and is 
only in svn as of today (13/5/2008)

In the unfortunate-but-common case of an unknown unit cell all is not lost. 
ImageD11 has a routine which will attempt to index an unknown unit cell 
from single crystal data. The script to use is called "index_unknown.py". 
As with most of the scripts, there is a list of options, which are more or 
less easy to set, depending on the problem. The options are::

 $ index_unknown.py --help
 Usage: index_unknown.py [options]

Options:
  -h, --help            show this help message and exit
  -g GVEFILENAME, --gve=GVEFILENAME
                        Filename for g-vectors
  -k NGRAINS, --ngrains=NGRAINS
                        number of grains to try to find
  -o OUTFILE, --output=OUTFILE
                        Name of ubi file to save grains in
  -v MIN_VEC2, --min_vec2=MIN_VEC2
                        Minimum axis length ^2, \AA^2 [1.5]
  -m N_TRY, --n_try=N_TRY
                        Number of vectors to test in finding lattice [all]
  -f FRACTION_INDEXED, --fraction_indexed=FRACTION_INDEXED
                        Fraction of peaks to be indexed
  -t TOL, --tol=TOL     tolerance in hkl error for indexing
  --fft                 Use fft to generate lattice vectors
  --score_fft           Score fft peaks using g-vectors
  --do_sort             Sorting the gvector by length before indexing [True]
  -n NP, --ngrid=NP     number of points in the fft grid [128]
  -r MR, --max_res=MR   Maximum resolution limit for fft (d-spacing) [1.0]
  -s NSIG, --nsig=NSIG  Number of sigma for patterson peaksearch threshold [5]

When you launch the script it will read in your g-vector file and then 
attempt to find a crystal lattice which accounts for the peaks in the 
g-vector file. It can attempt to generate the lattice either by combining 
g-vectors directly in reciprocal space, or alternatively by combining 
Patterson peaks from a fourier transform of the g-vector positions. The 
latter is interesting in the case of a dataset where the unit cell is large.

A more concrete explanation of the options follows::

 -g, --gve : The name of a file containing g-vectors
 -k, --ngrains : The number of grains you expect to find in the dataset. Only tested up to 3. It is really not intended for many grains just yet.
 -o --output : The ubi file which will receive the orientation matrices if and when they are found
 -v, --min_vec2 : Disturbingly non-obvious, and buggy still?. Should be the axis length when --fft is supplied or a g-vector error when using g-vectors. It is how close to zero a vector should be to be ignored when building a lattice. Feel free to edit the code to make this better, but make the testcases pass before committing to svn please, also edit here!
 -m, --n_try : When generating lattices from vectors you can use all possible combinations of choosing 3 vectors from all possible. That is often a large number. To avoid the problem we use only the first "n" with vectors sorted by length (gv) or peak height (patterson)
 -f, --fraction_indexed : Should be something like 1/k for now. TODO: make this a completeness or take account of k too?
 -t, --tol : The indexing error on g-vectors before you call them indexed (hkl units)
 -fft : A logical flag. By default the program tries to combine g-vectors into lattices. Adding --fft means it will use the fft
 --score_fft : Unlikely to be useful, but added for completeness (TODO: testcase missing). Scores how many peaks a trial unit cell indexes from the fft peaks instead of the g-vectors. May be much faster for larger numbers of peaks.
 --do_sort : Likely that you want it to be true. It decides whether to sort the g-vectors by length. Normally ImageD11 will have done this during the transform stage, but in certain cases it is needed to the n_try optimisation
 -n, --ngrid : Number of points in the FFT grid. 128 seems to be good, should be a multiple of 2. 256 takes significantly longer
 -r, --max_res : The maximum resolution (cut off) where peaks are put into the fft. Should be a d-spacing in angstrom. This determines the resolution of the fft.
 -s, --nsig : Threshold for peaksearching the patterson. Numbers like 10-30 seem to be useful.

Six synthetic testcases are supplied with the program in 
test/test_index_unknown/test_index_unknown.py. These generate pseudo-data 
g-vectors and then index them using the script. Firstly there are two 
tests of a single unknown via the gector and fft methods. Then there are 
two tests with two unknowns for the two methods. Finally three unknowns 
with the two methods. In the current SVN an issue remains that the program 
has a tendency to find "superlattices" in the case of several unknowns 
which somehow index more than one cell at the same time. This point bears 
further investigation, perhaps it is a problem with the choices of testcase.

As always, feedback and improvements are most welcome, please mail Jon 
directly if you don't get it to work first time.
