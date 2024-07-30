#!/usr/bin/env bash


# pip install dist/ImageD11-1.7.0-cp27-none-linux_x86_64.whl --user
# PYTHONPATH=$HOME/.local/lib/python2.7/site-packages
PYT=python
SRC=`pwd`
$PYT setup.py build_ext --inplace > bld.log 2> bld.err
PYTHONPATH=`pwd`
export PYTHONPATH=$PYTHONPATH
echo "Running tests from " $SRC " with PYTHONPATH: " $PYTHONPATH

cd $SRC/test && $PYT -c 'import ImageD11, sys; sys.stdout.write(ImageD11.__file__+"\n")'


cd $SRC/test
$PYT run_tests.py
cd $SRC

cd $SRC/test/demo
echo `pwd` test.py
$PYT test.py

cd $SRC/test/quantix/
echo `pwd` testfitfilt.py
$PYT testfitfilt.py 

cd $SRC/test
echo `pwd` test_peaksearch.py
$PYT test_peaksearch.py

cd $SRC/test
## takes really ages
#python2.5 test_peaksearch.py ALL

cd $SRC/test/ken_simul
$PYT testken.py


cd $SRC/test/
$PYT testcol.py


echo
echo "Just finished testing ImageD11 from" $PYT
echo "Using PYTHONPATH=" $PYTHONPATH
cd $SRC/test && $PYT -c 'import ImageD11, sys; sys.stdout.write(ImageD11.__version__+" "+ImageD11.__file__+"\n")'
cd $SRC
