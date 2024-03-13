from __future__ import print_function, division

inp = """
#File names and Format
#Directory to save output files
direc 'Al1000'
stem  'Al1000'
make_image 0
output '.gve' '.flt' '.par' '.ubi' 
# Structural parameters
structure_phase_0 'Al.cif'
# Crystal/grains parameters
no_grains 1000
# Total number of grains summed over all phases to be simulated. 
# Need to match the number of U_Grains_X key word
gen_size 1 -0.05 0.01 0.2
#Grain phase : If you want to let the PolyXSim appoint which grain 
# belongs to which phase the following keyword can be used.
gen_phase 1 1 0
# Grain orientation : random (1) or specify specific orientation matrices (0)
gen_U 0

# 1 random, box or cylinder
gen_pos 0 0 
# pos_grains_0 0.1 -0.1 0.05

#sample_xyz 1.0 1.0 1.0
gen_eps 1 0 0 0 0

# Instrumentation
# Detector pixel size [mm]
y_size 0.05 
z_size 0.05
#Detector size [pixels]
dety_size 2048
detz_size 2048
#Distance from sample to detector [mm]
distance 300
#Detector tilt 
tilt_x 0.005
tilt_y 0.01
tilt_z 0.008
#Detector orientation 
o11 1
o12 0
o21 0
o22 -1
#Noise
noise 0
#Reflection
intensity_const 1
lorentz_apply 1
beampol_apply 1
peakshape 0
#Instrumental
#Beam specs, Pt edge
wavelength 0.158154
beamflux 1e-12
beampol_factor 1
beampol_direct 0
#Beam center on detector [pixels]
dety_center 1022.1
detz_center 1028.3
#Omega scan range [degrees]
omega_start 0
omega_end 180
#Omega step size [degrees]
omega_step 0.25
#Omega rotation direction [degrees]
omega_sign 1
#Wedge angle of omega axis
wedge 0.02
"""

"U_grains_0 7.712806e-01 -6.337184e-01 5.939117e-02 6.130920e-01 7.146102e-01 -3.368241e-01 1.710101e-01 2.961981e-01 9.396926e-01"

import numpy as np, xfab.tools

# generate a 10x10x10 grid of grains.
# add something to this

x,y,z = np.mgrid[0:10,0:10,0:10]-5

dx = np.sin(np.arange(1000))/3 # +/- 1
dy = np.cos(np.arange(1000))/5
dz = np.sin(np.arange(1000))/7

t = np.array( (x.ravel()+dx, y.ravel()+dy, z.ravel()+dz ) )/10
#print t.shape
np.savetxt("t",t.T)

# orientations....
# t is in range [-0.5 -> 0.5] - make it rod also.

u = [xfab.tools.rod_to_u( v) for v in t.T] 

f=open("Al1000.inp","w")
f.write(inp)
for i,v in enumerate(t.T):
    f.write("pos_grains_%d %f %f %f\n"%(i,v[0],v[1],v[2]))
for i,v in enumerate(u):
    f.write("U_grains_%d "%(i))
    f.write(" ".join( repr(x) for x in v.ravel() ) )
    f.write('\n')
