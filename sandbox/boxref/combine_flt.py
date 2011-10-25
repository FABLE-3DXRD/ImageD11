

import columnfile, os, sys, glob, numpy

def sand_special_case(fname):
    """
    Create an ascii file holding fltname, x, y, z in microns
    one file per line
    """
    rootdir = "/data/id11/inhouse/jon/sand0910/newpeaksearches/peaksearches"
    stem = "sand_097_0910"
    scans = "a_","b_"
    xr = range(7)
    zr = range(1,2)
    fltnames = [ ("%s%s%d_%d__all.flt"%(stem, s, xp, zp),xp,zp)
                 for s in scans
                 for zp in zr
                 for xp in xr ]
    #for f,x,z in fltnames:
    #    print f, x, z, os.path.exists(os.path.join(rootdir, f))
    assert len(fltnames) == len(glob.glob(os.path.join(rootdir,"*_all.flt")))
    f = open(fname,"w")
    y = 0
    for flt, x,z in fltnames:
        f.write("%s %f %f %f\n"%(
            os.path.join(rootdir, flt),
            x * 1500 - 4500,
            y,
            z * 500 - 750
            ))
    f.close()
    # print open(fname).read()
    print "Hardwired sand special case, using default file with fltname samtx samty samtz"


class combine_flt(object):
    def __init__(self, fltfile, samtx=0, samty=0, samtz=0):
        print "reading", fltfile
        sys.stdout.flush()
        self.fltfile = columnfile.columnfile( fltfile )
        onearray = numpy.ones(self.fltfile.nrows, numpy.float32)
        self.fltfile.addcolumn( onearray*samtx, "samtx" )
        self.fltfile.addcolumn( onearray*samty, "samty" )
        self.fltfile.addcolumn( onearray*samtz, "samtz" )
        self.datalist = [ list(col) for col in self.fltfile.bigarray ]
        print len(self.datalist),len(self.datalist[0])
    def addfltfile(self, fltfile, samtx=0, samty=0, samtz=0):
        c = columnfile.columnfile(fltfile)
        onearray = numpy.ones(c.nrows, numpy.float32)
        c.addcolumn( onearray*samtx, "samtx" )
        c.addcolumn( onearray*samty, "samty" )
        c.addcolumn( onearray*samtz, "samtz" )
        for i,j in zip(c.titles, self.fltfile.titles):
            assert i == j
        for i in range(len(self.datalist)):
            self.datalist[i] += list(c.bigarray[i])
        print len(c.bigarray[0]), len(self.datalist[0])
    def write_hdf(self, name):
        self.datalist = numpy.array(self.datalist)
        self.fltfile.bigarray = self.datalist
        self.fltfile.ncols, self.fltfile.nrows = self.fltfile.bigarray.shape
        self.fltfile.set_attributes()
        columnfile.colfile_to_hdf(self.fltfile, name, "combinefltmergedpeaks")
        
        
if __name__ == "__main__":
    if len(sys.argv)==1:
        fname="sand_flt_files.txt"
        sand_special_case(fname)
    else:
        print "Reading fltname x y z from your file", fname
    o = None
    for line in open(fname).readlines():
        name = line.split()[0]
        x,y,z = [float(v) for v in line.split()[1:]]
        print name, x, y, z
        if o is None:
            o = combine_flt( name, x, y, z )
        else:
            o.addfltfile( name, x, y, z)
    o.write_hdf( "trial.hdf"  )
