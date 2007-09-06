#!/usr/bin/env sh

SRC="/home/wright/workspace/ImageD11-trunk"

cd $SRC


python2.5 setup.py build

export PYTHONPATH=$SRC/build/lib.linux-i686-2.5/


cd $SRC/test/demo
python2.5 test.py

cd $SRC/test/quantix/
python2.5 fitgrain.py g3.pars g3.ubi g3.flt new.pars


cd $SRC/test/gv_general
python2.5 test_gv_general.py

python2.5 -c 'import ImageD11; print ImageD11.__version__'



cd $SRC
python2.4 setup.py build

export PYTHONPATH=$SRC/build/lib.linux-i686-2.4/


cd $SRC/test/demo
python2.4 test.py

cd $SRC/test/quantix/
python2.4 fitgrain.py g3.pars g3.ubi g3.flt new.pars


cd $SRC/test/gv_general
python2.4 test_gv_general.py

python2.4 -c 'import ImageD11; print ImageD11.__version__'

