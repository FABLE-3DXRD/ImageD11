


from ImageD11 import cImageD11, sparseframe
import numpy as np, fabio, h5py, os


def grab2d( im, dark, nsigma=3., blocksize=256 ):
    """ returns the fabio image adding 
    .corr = corrected data
    .sparseframe = segmented spots
    """
    data = np.empty( im.size, np.float32 ) # corrected data
    msk =  np.empty( im.shape, 'b' )  # mask
    csk =  np.empty( im.shape, 'b' )  # cleaned mask
    # subtract the background
    shape0 = im.shape
    data.shape = data.size
    dark.shape = dark.size
    im.shape = im.size
    cImageD11.uint16_to_float_darksub( data, dark, im )
    # compute the mean and std of the data trimming at nsigma
    avg, std = cImageD11.array_mean_var_cut( data, nsigma )
    threshold = avg + nsigma * std
    # remove the readout drift
    data.shape = im.size // blocksize, blocksize
    cImageD11.frelon_lines( data, threshold )
    data.shape = shape0
    dark.shape = shape0
    im.shape = shape0
    # threshold and clean
    nnz = cImageD11.make_clean_mask( data, threshold, msk, csk )
    row = np.empty( nnz, np.uint16 )
    col = np.empty( nnz, np.uint16 )
    tmp = np.empty( shape0[0], 'i' )
    cImageD11.mask_to_coo( csk, row, col, tmp )
    header = { 'threshold' : threshold,
               }
    spar = sparseframe.from_data_mask( csk, data, header )
    return spar, data


def test():
    img = fabio.open( '../Au6_s0_003_a/Au6_s0_003_a0001.edf.gz' ).data
    drk = fabio.open("../pks/bkg.edf").data
    return grab2d( img, drk, nsigma=2.5, blocksize=128 )

def do_img(fname, drk, nsigma=2.5, blocksize=128 ):
    imo = fabio.open( fname )
    imo.spar, imo.corr = grab2d( imo.data, drk, nsigma=2.5, blocksize=128 )
    return imo
    

def grab_folder(folder):
    ROOT= "/data/id11/nanoscope/Commissioning/2017Feb/gold"
    drk = fabio.open("../pks/bkg.edf").data
    h5 = h5py.File(os.path.join(ROOT, "Jan2020", folder+".hdf"))
    gr = h5.require_group(folder)
    for j in range(900):
        edf = os.path.join(ROOT, folder, folder+"%04d.edf.gz"%(j)) 
        if os.path.exists( edf ):
            imo = do_img( edf, drk )
            g = gr.require_group( str(j) )
            imo.spar.to_hdf_group( g )
            print(edf,imo.spar.nnz)
            
if __name__=="__main__":
    import multiprocessing
    folders=[]
    for i in range(1,83):
        for d in 'ab':
            folders.append(os.path.join("Au6_s0_%03d_"%(i)+d))

    with multiprocessing.Pool(8) as p:
        p.map( grab_folder, folders )


