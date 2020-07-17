

import fabio, numpy as np
from ImageD11.labelimage import labelimage
from ImageD11.peaksearcher import peaksearch
from ImageD11.blobcorrector import correctorclass, perfect

lines = [ ("tiftest%04d.tif"%(i),i-6) for i in range(6,17) ]

# give somes minimal options
thresholds = [1, 10, 100]
corrector  = perfect() # no spatial disortion
dims = fabio.open(lines[0][0]).data.shape
label_ims  = { t : labelimage( shape=dims,
                               fileout="peaks_t%d.flt"%( t ),
                               sptfile="peaks_t%d.spt"%( t ),
                               spatial=corrector )
               for t in thresholds }

for filename, omega in lines:
    frame = fabio.open(filename)
    frame.header['Omega'] = omega
    frame.data = frame.data.astype( np.float32 )
    # corrections like dark/flat/normalise would be added here
    peaksearch( filename, frame, corrector, thresholds, label_ims )

for t in thresholds:
    label_ims[t].finalise()
    


