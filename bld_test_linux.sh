#!/usr/bin/env bash

python setup.py build bdist_wheel

# pip install dist/ImageD11-1.7.0-cp27-none-linux_x86_64.whl --user
# PYTHONPATH=$HOME/.local/lib/python2.7/site-packages
PYT=python
SRC=`pwd`
PYTHONPATH=$SRC/build/lib.linux-x86-64-2.7:$PYTHONPATH

echo "Running tests from " $SRC " with PYTHONPATH " $PYTHONPATH




cd $SRC/test/demo
echo `pwd` latred_new.py
$PYT latred_new.py
echo `pwd` test.py
$PYT test.py


cd $SRC/test
python run_tests.py



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
$PYT idx.py

cd $SRC/test/
$PYT testcol.py


cd $SRC


echo
echo "Just finished testing ImageD11 from" $PYT
echo "Using PYTHONPATH=" $PYTHONPATH
$PYT -c "import ImageD11; print ImageD11.__version__"
