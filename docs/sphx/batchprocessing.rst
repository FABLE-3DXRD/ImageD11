===============
BatchProcessing
===============
Imagine you have collected a mass of data. You want to peaksearch it, 
using several machines, and then index all the scans. Here is one 
approach. 
These are notes from a worked example. You might need to either try to 
understand the code/method or wait for a more user friendly implementation


Working out what the scans were
===============================

This is an example script which will go through the headers getting 
image numbers and header information (svn version only for now)::

 more getheader.py
 from ImageD11 import opendata
 import glob, sys
 fl = glob.glob(sys.argv[1])   # list of files to process
 output = []
 for name in fl:
     if name.find("edf") > -1:   # only process edf files
         h = opendata.edfheader(name)
         num = opendata.getnum(name)  //**# this needs svn ImageD11.opendata**//
         for n in sys.argv[2:]:
             try:
                 name += " " + str(n) + " " + h[n] + " "
             except:
                 name += " " + str(n) + " " + " " + " "
         output.append((num,name))
 output.sort()
 for o in output:
     print o[0],o[1]

	 
Usage would be::

 python getheader.py "mydatafiles*.edf" samz Omega > headerinfo.txt

The quotes around the first command line argument are critical to stop the shell expanding the "*". 
It will create a file with the image number, filename and header information you want. We will put the Omega entry as the last column and last argument on the command line. We now work out what the individual sweepscans were from this header information file. It is sorted by image number and we assume a scan corresponds to an increasing sequence of Omega angles::

 % more getscans.py
 import sys
 njobs = 4   # number of shell scripts to create
 outfiles = [open("fable%d.sh"%(i+1),"w") for i in range(njobs)]
 last = -1e9
 first = True
 job = 0
 for line in open(sys.argv[1],"r").readlines():
     items = line.split()
     num = int(items[0])    # image number
     om = float(items[-1])  # Omega angle
     if first:
         f = num
         first = False
     if om > last:
         last = om
     else:
         # make a job
         stem = "/data/opid11/external/uu123/mysample/" # your directory here
         # build a peaksearch.py command line
         com = "peaksearch.py "
         com += " -n %s "%(stem+"mysample")
         com += " -f %d "%(f)
         com += " -l %d "%(num)
         # we made bkg.edf via ImageD11/scripts/bgmaker.py
         com += " -d bkg.edf "
         com += " -D 20 "
         com += " -t 70 -t 520 -t 5020 "
         com += " -s /data/opid11/inhouse/Frelon4M/frelon4m.spline "
         # output to go in ImageD11/peaksearch subdirectory
         com += " -o %sImageD11/peaksearch/scan_%d.out "%(stem,job)
         com += " >  %sImageD11/peaksearch/scan_%d.log "%(stem,job)
         com += "\n"
         outfiles[job%njobs].write("echo %d\n"%(job))
         outfiles[job%njobs].write(com)
         job += 1
         first = True
         last = -1e9
		 
Assuming you have made a directory ImageD11 in with your sample to pick up the output should now have 4 shell scripts. Run each one on a different machine, the peaksearch output should end up in the ImageD11 directory. 

Work out how to index one scan
==============================

Take the output from the peaksearch into the gui (probably the merge_txxx) file. Go through the transformation menu until you are happy you have nice parameters. Save the parameter file. Save the g-vectors and go through the indexing menu, finding grains. Once you are happy go to the Help menu and select history. Copy and paste this into a file.

Edit the history file to take a command line argument
=====================================================

Now edit the history file to replace the input and output file names to be a variable from the command line: eg::

 % more history.py
 # Create objects to manipulate - they hold your data
 #
 from ImageD11 import peakmerge, indexing, transformer
 mypeakmerger = peakmerge.peakmerger()
 mytransformer = transformer.transformer()
 myindexer = indexing.indexer()
 # 
 # Your work starts here:
 #
 **//import sys  # added//**
 //**infile = sys.argv[1] # added**//
 mytransformer.loadfiltered( //**infile**// )
 mytransformer.loadfileparameters(  '/data/opid11/external/uu123/mysample/ImageD11/good.pars' )
 mytransformer.compute_tth_eta( )
 mytransformer.addcellpeaks( )
 mytransformer.savegv( //**infile+'.gve'  ) # add .gve extension**//  

 myindexer.parameterobj.set_parameters(
  {'ds_tol': '0.01',
   'minpks': '50',
   'uniqueness': '0.5',
   'hkl_tol': '0.05',
   'eta_range': '10.0',
   'ring_1': '3',
   'wavelength': '0.29173',
   'ring_2': '3',
   'cosine_tol': '0.005'} )
 myindexer.updateparameters( )
 myindexer.loadpars( )
 myindexer.readgvfile(  //**infile+'.gve'**// )
 myindexer.assigntorings( )
 myindexer.find( )
 myindexer.scorethem( )
 myindexer.assigntorings( )
 myindexer.saveubis(   //**infile + '.ubi'**// )
 //**# slow process, skip to just get ubi files**//
 //**# myindexer.saveindexing( infile + '.idx')**//

 Now run this on all you peaksearch outputs:

 bash % find . -name "*_merge_t500" -exec history.py {} \;

 And off it goes to make the ubi files for you

Indexing a peaksearch output and filtering peaks at the same time::

 # Example indexing script. Typically will be edited to match
 # specific experiments and sample (eg: which rings to use for indexing,
 # how many peaks make a grain, how to filter the peaks to remove "bad"
 # peaks coming from the peaksearch)
 #
 # call as idx.py [pars] [flt] [gve]  [ubi]
 #                input  input output output
 #
 # Create objects to manipulate - they hold your data
 #
 import sys
 from ImageD11 import peakmerge, indexing, transformer
 mypeakmerger = peakmerge.peakmerger()
 mytransformer = transformer.transformer()
 myindexer = indexing.indexer()
 #
 # Your work starts here:
 #
 mytransformer.loadfiltered( sys.argv[2] )
 mytransformer.loadfileparameters( sys.argv[1] )
 mytransformer.compute_tth_eta( )
 # Cut the high angle rings, so only tth less than 14
 # ... could be an option
 mytransformer.colfile.filter( mytransformer.colfile.tth < 14 )
 # ... you can filter on any column here
 # eg   ... mytransformer.colfile.Number_of_pixels > 2
 #      ... etc
 mytransformer.compute_tth_eta( )
 mytransformer.addcellpeaks( )
 mytransformer.computegv( )
 mytransformer.savegv( sys.argv[3] )
 
 
 myindexer.readgvfile( sys.argv[3] )
 myindexer.updateparameters( )
 # Mostly defaults - wavelength is common to pars above
 myindexer.parameterobj.set_parameters(  {
      'ds_tol': '0.02',
      'minpks': '30',
      'max_grains': '100',
      'uniqueness': '0.5',
      'hkl_tol': '0.07',
      'eta_range': '0.0',
      'ring_1': '1',
      'ring_2': '1',
      'cosine_tol': '0.01'} )
 myindexer.loadpars( )
 myindexer.assigntorings( )
 myindexer.find( )
 myindexer.scorethem( )
 
 # loop over minpks reducing from 40 to 20 - gets the best grains
 # with 40 peaks first
 for minpk in [ 40, 30, 20]:
      # loop over pairs of rings to use for generating orientations
      for p1, p2 in [  (0, 1), (2, 1), (3, 1) , (4, 1) ,
                       (0, 0), (2, 0), (3, 0) , (4, 0) ]:
 
          myindexer.updateparameters( )
          myindexer.parameterobj.set_parameters(  {
              'ring_1': p1,
              'ring_2': p2,
              'minpks': minpk,
              } )
          myindexer.loadpars( )
          myindexer.assigntorings( )
          myindexer.find( )
          myindexer.scorethem( )
 
 myindexer.saveubis( sys.argv[4] ) 
