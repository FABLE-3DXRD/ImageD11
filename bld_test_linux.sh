#!/usr/bin/env bash



SRC=`dirname "$0"`
PYTHONPATH=$SRC/build/lib.linux-i686-2.5/

echo "Running tests from " $SRC " with PYTHONPATH " $PYTHONPATH

# Hmm - needs a switch case

case `hostname`  in
    lapwright) PYT=python2.5 ;;
    *)  PYT=fable.python ;;
esac

cd $SRC
$PYT setup.py build

echo test_put_incr.py
cd $SRC/test
$PYT test_put_incr.py


cd $SRC/test/demo
echo `pwd` latred_new.py
$PYT latred_new.py
echo `pwd` test.py
$PYT test.py

cd $SRC/test/test_index_unknown
echo `pwd` test_index_unknown.py
$PYT test_index_unknown.py

cd $SRC/test/quantix/
echo `pwd` testfitfilt.py
$PYT testfitfilt.py 

cd $SRC/test
echo `pwd` test_peaksearch.py
$PYT test_peaksearch.py

cd $SRC/test
## takes really ages
#python2.5 test_peaksearch.py ALL


cd $SRC/test/gv_general
$PYT test_gv_general.py

cd $SRC/test/testconnectedpixels
$PYT testconnectedpixels.py
 
cd $SRC/test/testlabelimage
$PYT testlabelimage.py

cd $SRC/test/ken_simul
$PYT idx.py

cd $SRC/test/
$PYT test_ubito.py
$PYT testcolumnfile.py
$PYT testscale.py
$PYT testcol.py


cd $SRC


echo
echo "Just finished testing ImageD11 from" $PYT
echo "Using PYTHONPATH=" $PYTHONPATH
$PYT -c "import ImageD11; print ImageD11.__version__"
