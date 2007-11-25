

if [ $HOSTTYPE = 'i686' ]; then
export PYTHONPATH=../build/lib.linux-i686-2.5
export LD_LIBRARY_PATH=/sware/exp/fable/standalone/suse82/lib
/sware/exp/fable/standalone/suse82/bin/python2.5  << EOF
 

from ImageD11 import blobcorrector
b = blobcorrector.correctorclass("spatial2k.spline")
print b.correct(10,20)
EOF

fi





if [ $HOSTTYPE = 'x86_64' ]; then
export PYTHONPATH=../build/lib.linux-x86_64-2.5
export LD_LIBRARY_PATH=/sware/exp/fable/standalone/redhate4-a64/lib

/sware/exp/fable/standalone/redhate4-a64/bin/python2.5  << EOF
 

from ImageD11 import blobcorrector
b = blobcorrector.correctorclass("spatial2k.spline")
print b.correct(10,20)
EOF

fi
