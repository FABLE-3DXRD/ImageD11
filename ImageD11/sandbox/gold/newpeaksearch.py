
import pylab as pl, numpy as np, h5py, time, glob
import multiprocessing
from ImageD11 import sparseframe, cImageD11


def newpks(hf):
    """ do a peak search in a sparse frame series """
    blbs = []
    with h5py.File(hf,"r") as h:
        stem = list(h)[0]
        grp = h[stem]
        for g in range(900):
            f = sparseframe.from_hdf_group(grp[str(g)])
            sparseframe.sparse_localmax(f)
            blbs.append( sparseframe.sparse_moments(f, "intensity",
                                                    "localmax") )
    return stem, blbs


def main():
    ################### scan specific motorpos ########
    dtyomepos = {}
    oa = np.arange(0.,180.,0.2)
    ob = np.arange(180.,0.,-0.2)
    assert len(oa)==len(ob)

    with open( "dty.pos", "r") as f:
        for line in f.readlines()[1:]:
            scan = line.split("/")[0]
            dty = float(line.split()[1])
            if scan[-1] == 'a':
                dtyomepos[scan]=dty, oa
            if     scan[-1] == 'b':
                dtyomepos[scan]=dty, ob
        


    hdfs = [ "Au6_s0_%03d_%s.hdf"%(i,ab)
             for  i in range(1,81) 
             for ab in 'ab'] 

    pks = {}
    start = time.time()
    with multiprocessing.Pool(multiprocessing.cpu_count()-1) as p:
        for ret, h in zip(p.imap_unordered(newpks, hdfs), hdfs):
            s,b = ret
            d,o = dtyomepos[s]
            pks[d] = o, b
            end = time.time()
            print(h,d,o[0],len(b),end-start)




if __name__=="__main__":
    main()
