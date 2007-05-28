


from ImageD11 import opendata
import glob

def getnumerics(d):
    ret = {}
    for k,v in d.items():
        try:
            ret[k] = float(v)
        except ValueError:
            pass
    return ret

import sys
class progress:
    def __init__(self):
        self.whirl = [ "-\r", "\\\r", "|\r", "/\r" ]
        self.i = 0
        self.l = len(self.whirl)
    def update(self):
        self.i += 1
        self.i %= self.l
        sys.stdout.write(self.whirl[self.i])
        sys.stdout.flush()
        
    

class headerfollower:
    
    # Header items we choose to ignore changes in for printing
    
    ignore = ['time_of_day','stream700','pico1','Ring Current','samome',
              'run']
    interesting = [ 'time' ,  'dir', 'prefix', 'suffix']

    
    def __init__(self, filename):
        """
        Initialise on a filename
        """
        self.progress = progress()
        self.header = opendata.edfheader(filename)
        self.numerics = getnumerics(self.header)
        self.headeritems = set(self.header)
        self.filename = filename
        self.report = {}
        self.reportnames = []
        for k in self.interesting:
            self.addreport(filename, k, self.header[k])
        for h,v in self.numerics.items():
            self.addreport(filename, h, v)
        
    def addreport(self, name, item, val):
        """
        Store up the things to report in a dictionary
        Holds a list of names, then item val pairs to report
        """
        if self.report.has_key(name):
            self.report[name].append( (item, val) )
        else:
            self.report[name] = [(item, val)]
            self.reportnames.append(name)

    def printreport(self):
        """
        pretty print the report information
        """
        print
        for name in self.reportnames:
            print name
            sr = self.report[name]
            sr.sort()
            for k, v in sr:
                print "\t",k,v
            print
                                       
    def nextfile(self,filename):
        """
        process the next file
        """
        self.progress.update()
        h = opendata.edfheader(filename)
        hks = set(h.keys())
        changes = hks.symmetric_difference(self.headeritems)
        if len(changes) != 0:
            raise Exception("Filename %s introduces different header keys"%
                            (filename))
        n = getnumerics(h)
        for k in n.keys():
            if k == "Omega":
                if n[k] < self.numerics[k]:
                    self.addreport(self.filename,k,self.header[k])
                    self.addreport(filename, k, h[k])
                self.numerics[k] = n[k]
                self.header[k] = h[k]
                continue
            if k in self.ignore:
                continue
            if n[k] != self.numerics[k]: # number changed
                # report
                self.addreport(self.filename,k,self.header[k])
                self.addreport(filename, k, h[k])
                self.numerics[k] = n[k]
                self.header[k] = h[k]
        self.filename = filename
        
            

filelist = glob.glob("*.edf")

filelist = [ (opendata.getnum(name), name) for name in filelist ]

filelist.sort()

obj = headerfollower(filelist[0][1])

for n,f in filelist[1:]:
    obj.nextfile(f)

obj.printreport()




