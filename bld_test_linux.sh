#!/usr/bin/env bash



# Hmm - needs a switch case

if [ `hostname` = "lapwright" ]; then
    SRC="/home/wright/workspace/ImageD11-trunk"
    PYT=python2.5
else
    if [ "$HOSTTYPE" = "x86_64" ]; then
	SRC="/users/wright/fable/svn/ImageD11/trunk"
	export LD_LIBRARY_PATH=/sware/exp/fable/standalone/redhate4-a64/lib
	PYT=/sware/exp/fable/standalone/redhate4-a64/bin/python
    fi
fi




cd $SRC
$PYT setup.py build
export PYTHONPATH=$SRC/build/lib.linux-i686-2.5/
PYTHONPATH=$SRC/build/lib.linux-i686-2.5/


cd $SRC/test
$PYT test_put_incr.py


cd $SRC/test/demo
$PYT latred_new.py
$PYT test.py


cd $SRC/test/quantix/
$PYT fitgrain.py g3.pars g3.ubi g3.flt new.pars

cd $SRC/test
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



echo
echo "Just finished testing ImageD11 from" $PYT
echo "Using PYTHONPATH=" $PYTHONPATH
$PYT -c "import ImageD11; print ImageD11.__version__"
