from __future__ import print_function

import os
import ImageD11.sinograms.dataset

    
#    dataroot = sys.argv[1]
#    analysisroot = sys.argv[2]
#    sample = sys.argv[3]
#    dset = sys.argv[4]
#    destination = sys.argv[5]


# need to locate some 'live' datasets for testing.

def testcase(dataroot = '/data/id11/nanoscope/ihmi1452/id11',
             sample = 'WAucross',
             dset = 'H0_',
             destination = '/tmp/case.h5',
             analysisroot = 'should_not_exist',
             scans= None
            ):
    dsname = "_".join( ( sample, dset ) )
    mfile = os.path.join( dataroot, sample, dsname, dsname+'.h5' )
    if not os.path.exists( mfile ):
        print( "Missing files", mfile )
        return # False
    try:
        result = ImageD11.sinograms.dataset.check( dataroot, 
                                                   analysisroot,
                                                   sample,
                                                   dset,
                                                   destination,
                                                   scans)
    except Exception as e:
        print("Error occured",e )
        raise
    finally:
        print("Remove",destination)
        os.unlink( destination )
    return # result




testcase(
    dataroot = '/data/id11/nanoscope/ihmi1452/id11',
    sample = 'WAucross',
    dset = 'H0_',
)

testcase(
    dataroot = '/data/id11/nanoscope/ihma195/id11',
    sample = 'goldBall',
    dset = 'r2_Z071',
    scans = [ '4.1',] 
)
