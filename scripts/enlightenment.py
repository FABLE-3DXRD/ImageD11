#!/usr/bin/env python

from __future__ import print_function



# python script to process the box scan style data
# Originally run as user in603 from
# /data/proposal/id11 via command "fable.python process_grid.py"

# Updated to be more general...




HOME = "/data/opid11/inhouse/gundlach"
JOBS = "resistor_map*"
SKIP = ["resistor_map_t0_000_"]


import sys, os, time, glob

class job:
    """ Represents a command line thing to do """
    def __init__(self, wk_dir, name, kwds):
        """
        wk_dir = working directory to use for running the command
        name   = a unique name for labelling a tee'd log file
        """
        self.wk_dir   = wk_dir
        self.name     = name
        self.kwds     = kwds
    def chk(self):
        """
        check if job is already run or not by looking for the log file
        """
        ret = os.path.exists(os.path.join(self.wk_dir, self.name+".log"))
        print("check",self.wk_dir, self.name+".log",ret)
        return ret
    def run(self):
        """ run the job as a command line """
        try:
            print(self.name, self.wk_dir)
            os.chdir(self.wk_dir)
            c = self.cmd() + " | tee %s.log"%(self.name)
            sys.stdout.write("running %s\n%s\n"%(self.name, c))
            ret = os.system(c)
            if ret != 0:
                print("Failed command")
                sys.exit()
            time.sleep(1)
        except:
            print(self.name, self.wk_dir, self.kwds)
            raise
    def cmd(self):
        """ this is the method to override"""
        raise Exception("Override me")


### classes for different jobs to do follow
    
class mkdir(job):
    """ Make a working directory for quantix / frelon """
    def cmd(self):
        return "mkdir %s"%(self.kwds['dirname'])

class bgmaker(job):
    """ background """
    def cmd(self):
        return "bgmaker.py -n %s -f %s -l %s -s %s -o %s"%(
            self.kwds['stem']  ,
            self.kwds['first'] ,
            self.kwds['last']  ,
            self.kwds['step']  ,
            self.kwds['bkg'] )

class peaksearch(job):
    """ peaksearch """
    def cmd(self):
        s = "peaksearch.py -n %s -f %s -l %s -d %s -o %s -p %s -s %s "%(
            self.kwds['stem']  ,
            self.kwds['first'] ,
            self.kwds['last']  ,          
            self.kwds['bkg']   ,
            self.kwds['out']   ,
            self.kwds['perfect'] ,
            self.kwds['spline'] )
        for t in self.kwds['thresholds']:
            s = s + " -t %s "%(t)
        return s    
        
class cat(job):
    """ join peaksearches together """
    def cmd(self):
        for f1,f2,res in zip(self.kwds['f1l'],
                             self.kwds['f2l'],
                             self.kwds['res']):
            os.system("cat %s %s > %s"% (f1, f2, res))
        return "echo catted"

class idx(job):
    """ index some grains """
    def cmd(self):
        return "fable.python %s/idx.py %s %s %s %s"%(HOME,
                                                     self.kwds['pars'],
                                                     self.kwds['flt'],
                                                     self.kwds['gve'],
                                                     self.kwds['ubi'])

class makemap(job):
    """ make / fit a grain map """
    def cmd(self):
        return "makemap.py -f %s -F %s -u %s -p %s -U %s -t %s -s cubic"%(
            self.kwds['flt'],
            self.kwds['not'],
            self.kwds['ubi'],
            self.kwds['par'],
            self.kwds['map'],
            self.kwds['tol'])
            
class collectmaps(job):
    def cmd(self):
        return 'cp %s %s'% (self.kwds['map'],
                            os.path.join(HOME,
                                         self.kwds['mapdir'],
                                         self.kwds['newname']))


            
def makejobs( direc, home ):
    """
    Generate the full list of jobs to do for each sample
    """
    qd = os.path.join(home, direc, "ImageD11_quantix")
    fd = os.path.join(home, direc, "ImageD11_frelon")
    print("direc",direc)
    print("qd", qd)
    print("fd",fd)
    joblist = [
        mkdir( os.path.join(home, direc) ,
               "mkdir_quantix", {'dirname' : qd } ) ,
        mkdir( os.path.join(home, direc) ,
               "mkdir_frelon" , {'dirname' : fd } ) ,

        bgmaker( qd ,  "bkg_quantix_a",{ 'stem'  : os.path.join("..",direc+"_quantix_") ,
                                         'first' : 0   ,
                                         'last'  : 179 ,
                                         'step'  : 3   ,
                                         'bkg'   : "bkg_qa.edf" }),
        
        peaksearch(qd, "pksh_a",      { 'stem'  : os.path.join("..",direc+"_quantix_") ,
                                        'first' : 0   , 
                                        'last'  : 179 ,
                                        'bkg'   : "bkg_qa.edf",
                                        'out'   : "peaks_qa.spt",
                                        'perfect':'Y',
                                        'spline': None,
                                        'thresholds' : [50, 100, 500, 3000] }),

        bgmaker( qd ,  "bkg_quantix_b",{'stem'  : os.path.join("..",direc+"_quantix_") ,
                                        'first' : 180 ,
                                        'last'  : 359 ,
                                        'step'  : 3   ,
                                        'bkg'   : "bkg_qb.edf"}),

        peaksearch(qd, "pksh_b",       {'stem'  : os.path.join("..",direc+"_quantix_") ,
                                        'first' : 180 ,
                                        'last'  : 359 ,
                                        'bkg'   : "bkg_qb.edf",
                                        'out'   : "peaks_qb.spt",
                                        'perfect':'Y',
                                        'spline': None,
                                        'thresholds' : [50, 100, 500, 3000]} ),
        

        bgmaker( fd ,  "bkg_frelon_a", {'stem'  : os.path.join("..",direc+"_frelon4m_"),
                                        'first' : 0   ,
                                        'last'  : 179 ,
                                        'step'  : 3   ,
                                        'bkg'   : "bkg_fa.edf"}),

        peaksearch(fd, "pksh_a",       {'stem'  : os.path.join("..",direc+"_frelon4m_"),
                                        'first' : 0   ,
                                        'last'  : 179 ,
                                        'bkg'   : "bkg_fa.edf",
                                        'out'   : "peaks_fa.spt",
                                        'perfect':'N',
                                        'spline': "/data/opid11/inhouse/Frelon4M/frelon4m.spline",
                                        'thresholds' : [200, 1000, 10000] }),

        bgmaker( fd ,  "bkg_frelon_b", {'stem'  : os.path.join("..",direc+"_frelon4m_"),
                                        'first' : 180 ,
                                        'last'  : 359 ,
                                        'step'  : 3   ,
                                        'bkg'   : "bkg_fb.edf"}),

        peaksearch(fd, "pksh_b",       {'stem'  : os.path.join("..",direc+"_frelon4m_"),
                                        'first' : 180 ,
                                        'last'  : 359 ,
                                        'bkg'   : "bkg_fb.edf",
                                        'out'   : "peaks_fb.spt",
                                        'perfect':'N',
                                        'spline': "/data/opid11/inhouse/Frelon4M/frelon4m.spline",
                                        'thresholds' : [200, 1000, 10000] }),

        cat( qd, "cat_2scans" , { 'f1l' : ["peaks_qa_t50.flt"   ,
                                           "peaks_qa_t100.flt"  ,
                                           "peaks_qa_t500.flt"  ,
                                           "peaks_qa_t3000.flt"],
                                  'f2l' : ["peaks_qb_t50.flt"   ,
                                           "peaks_qb_t100.flt"  ,
                                           "peaks_qb_t500.flt"  ,
                                           "peaks_qb_t3000.flt"],
                                  'res' : ["t50.flt"  ,
                                           "t100.flt" ,
                                           "t1500.flt",
                                           "t3000.flt"  ] } ),

        cat( fd, "cat_2scans" , { 'f1l' : ["peaks_fa_t200.flt",
                                           "peaks_fa_t1000.flt",
                                           "peaks_fa_t10000.flt"],
                                  'f2l' : ["peaks_fb_t200.flt",
                                           "peaks_fb_t1000.flt",
                                           "peaks_fb_t10000.flt"],
                                  'res' : ["t200.flt",
                                           "t1000.flt",
                                           "t10000.flt"] } ),

        idx( fd, "idx_10k" ,  { 'pars' : os.path.join( HOME, "frelon4m_resistor.par"),
                                'flt'  : 't10000.flt',
                                'gve'  : 't10000.gve',
                                'ubi'  : 't10000.ubi' }),


        makemap( fd, "map_10k1", { 'flt' : "t10000.flt",
                                  'not' : "not10000.flt",
                                  'ubi' : "t10000.ubi",
                                  'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                  'map' : "map10000.ubi",
                                  'tol' : 0.08 } ),

        makemap( fd, "map_10k2", { 'flt' : "t10000.flt",
                                  'not' : "not10000.flt",
                                  'ubi' : "map10000.ubi",  # <--- change here to recycle
                                  'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                  'map' : "map10000.ubi",
                                  'tol' : 0.04 } ) ,        # <---

        makemap( fd, "map_10k3", { 'flt' : "t10000.flt",
                                  'not' : "not10000.flt",
                                  'ubi' : "map10000.ubi",  # <--- change here to recycle
                                  'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                  'map' : "map10000.ubi",
                                  'tol' : 0.025 } ),         # <---
        
        # use the 10k ubis to filter the 1k
        makemap( fd, "map_1k",  { 'flt' : "t1000.flt",
                                  'not' : "not1000.flt",
                                  'ubi' : "map10000.ubi",  # <--- change here to recycle
                                  'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                  'map' : "map1000.ubi",
                                  'tol' : 0.025 } ),         # <---
        
        idx( fd, "idx_1k_not" ,  { 'pars' : os.path.join( HOME, "frelon4m_resistor.par"),
                                'flt'  : 'not1000.flt',
                                'gve'  : 'not1000.gve',
                                'ubi'  : 'not1000.ubi' }),

        # Fit alone
        makemap( fd, "map_1k_not1",  { 'flt' : "not1000.flt",
                                     'not' : "/dev/null",
                                     'ubi' : "not1000.ubi",
                                     'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                     'map' : "not1000.ubi",
                                     'tol' : 0.08 } ),         # <---
        makemap( fd, "map_1k_not2",  { 'flt' : "not1000.flt",
                                     'not' : "/dev/null",
                                     'ubi' : "not1000.ubi", 
                                     'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                     'map' : "not1000.ubi",
                                     'tol' : 0.04 } ),         # <---
        
        makemap( fd, "map_1k_not3",  { 'flt' : "not1000.flt",
                                     'not' : "/dev/null",
                                     'ubi' : "not1000.ubi",
                                     'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                     'map' : "not1000.ubi",
                                     'tol' : 0.025 } ),         # <---

        cat( fd, "cat_10_1_ubi" , { 'f1l' : ['map1000.ubi' ],
                                    'f2l' : ['not1000.ubi' ],
                                    'res' : ["final_1000.ubi"] } ),
        
        makemap( fd, "final_1k",   { 'flt' : "t1000.flt",
                                     'not' : "not1000.flt",
                                     'ubi' : "final_1000.ubi", 
                                     'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                     'map' : "final_1000.ubi",
                                     'tol' : 0.025 } ),
        # repeat...
        idx( fd, "idx_1k_not2" ,  { 'pars' : os.path.join( HOME, "frelon4m_resistor.par"),
                                'flt'  : 'not1000.flt',
                                'gve'  : 'not1000.gve',
                                'ubi'  : 'not1000.ubi' }),


        makemap( fd, "map_1k_not1b",  { 'flt' : "not1000.flt",
                                       'not' : "/dev/null",
                                       'ubi' : "not1000.ubi",
                                       'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                       'map' : "not1000.ubi",
                                       'tol' : 0.08 } ),         # <---
        makemap( fd, "map_1k_not2b",  { 'flt' : "not1000.flt",
                                       'not' : "/dev/null",
                                       'ubi' : "not1000.ubi", 
                                       'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                       'map' : "not1000.ubi",
                                       'tol' : 0.04 } ),         # <---
        
        makemap( fd, "map_1k_not3b",  { 'flt' : "not1000.flt",
                                        'not' : "/dev/null",
                                        'ubi' : "not1000.ubi",
                                        'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                        'map' : "not1000.ubi",
                                        'tol' : 0.025 } ),         # <---

        cat( fd, "cat_10_1_ubi" , { 'f1l' : ['final_1000.ubi' ],
                                    'f2l' : ['not1000.ubi' ],
                                    'res' : ["map_1000.ubi"] } ),
        
        makemap( fd, "final_1k",   { 'flt' : "t1000.flt",
                                     'not' : "not1000.flt",
                                     'ubi' : "map_1000.ubi", 
                                     'par' : os.path.join( HOME, "frelon4m_resistor.par"),
                                     'map' : "final_1000.ubi",
                                     'tol' : 0.025 } ),
        
        # kill bad   ??

        
        makemap( qd, "map_q1", { 'flt' : "t100.flt",
                                 'not' : "/dev/null",
                                 'ubi' : os.path.join("..","ImageD11_frelon","final_1000.ubi"),
                                 'par' : os.path.join( HOME, "quantix_resistor.par"),
                                 'map' : 'qmap100.ubi',
                                 'tol' : 0.15 } ),
        
        makemap( qd, "map_q2", { 'flt' : "t100.flt",
                                 'not' : "/dev/null",
                                 'ubi' : 'qmap100.ubi',
                                 'par' : os.path.join( HOME, "quantix_resistor.par"),
                                 'map' : 'qmap100.ubi',
                                 'tol' : 0.10 } ),

        makemap( qd, "map_q3", { 'flt' : "t100.flt",
                                 'not' : "/dev/null",
                                 'ubi' : 'qmap100.ubi',
                                 'par' : os.path.join( HOME, "quantix_resistor.par"),
                                 'map' : 'qmap100.ubi',
                                 'tol' : 0.05 } ),

        # Score/ check them
        # ... other thresholds ??

        collectmaps( qd , 'collect_resistor_maps',
                     {  'map':'qmap100.ubi',
                        'mapdir' : 'resistor_qmaps',
                        'newname': os.path.join( HOME, 'resistor_qmaps',
                                                 direc + 'qmap100.ubi') })

                                              

        ]
    return joblist

       



# wait on new directories to appear
# loop over directories doing all jobs for each one
if 1:
    # Return to start dir 
    os.chdir(HOME)

    dirlist = glob.glob(JOBS)
    dirlist.sort()

    # for each directory
    for d in dirlist:
        os.chdir(HOME)


        if not os.path.isdir(d) or d in SKIP:
            print("skipping",d)
            continue

        print("Entering",d)        

        
        jl = makejobs( d , HOME )
        for j in jl:
            os.chdir( os.path.join(HOME,d) )
            print(j.name, j.wk_dir, j.kwds)
            if not j.chk():
                print("running")
                j.run()
            else:
                print("done")
        
        
        

        
        
    
    
