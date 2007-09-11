



#!/usr/bin/env sh

SRC="/home/wright/workspace/ImageD11-trunk"

cd $SRC


python2.5 setup.py build

export PYTHONPATH=$SRC/build/lib.linux-i686-2.5/


#cd $SRC/test/demo
#python2.5 test.py

#cd $SRC/test/quantix/
#python2.5 fitgrain.py g3.pars g3.ubi g3.flt new.pars

#cd $SRC/test
#python2.5 test_peaksearch.py

cd $SRC/test/gv_general
python2.5 test_gv_general.py

cd $SRC/test/testconnectedpixels
python2.5 testconnectedpixels.py

cd $SRC/test/testlabelimage
python2.5 testlabelimage.py


python2.5 -c 'import ImageD11; print ImageD11.__version__'
