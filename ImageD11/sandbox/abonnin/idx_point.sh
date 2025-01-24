#!/bin/bash



cd /users/abonnin/lien_ma1063/diffraction/IRIS4_rect5_1/jon
source /sware/exp/fable/bin/fable_new.bash

python26 select_spatial_point.py mainphase_fixed_omega.flt  $1 $2 14.4 1.0  $3.flt
python26 makegv.py $3.flt mainphase.par $3.gve
cp template.ini $3.ini
echo filespecs $3.gve $3.log >> $3.ini
/sware/exp/fable/bin/GrainSpotter.0.90 $3.ini

python26 make_point_sino.py $3.log $3.flt > $3.sinfit

cp $3.sinfit gridsearchoutput
cp $3.gff gridsearchoutput

rm $3.flt $3.ini $3.gve $3.log $3.ubi $3.gff $3.sinfit



