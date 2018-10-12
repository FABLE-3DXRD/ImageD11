
from __future__ import print_function
"""
For a peak search interface
  Try to navigate the output of 'ls' in folders with lots of files
  Compress the list of names when there are a lot of numbers
"""


import re, os, time, sys
try:
    from os import scandir
except ImportError:
    class wrap:
        def __init__(self, name ):
            self.name = name
        def is_dir(self):
            return os.path.isdir( self.name )
    def scandir( folder ):
        for name in os.listdir( folder ):
            yield wrap( name )


            
def list_a_folder( folder ):
    """
    Lists files in the folder
      Groups into:
        directories
        names of type (stem)(num)(extn)
        things it doesn't understand (no number)
    """
#    NotGreedyAnything:Digit:~Digit:End
    reg = re.compile("(.*?)(\d*)(\D*)$")
    direcs = []
    names = []
    files  = {} # keyed by extension, then stem, then num sorted
    for f in scandir( folder ):
        if f.is_dir():
            direcs.append( f.name )
            continue
        items = reg.match( f.name )
        if items:
            stem, num, extn = items.groups()
            if len(num)==0:
                names.append( f.name )
                continue
            if extn in files:
                if stem in files[extn]:
                    files[extn][stem].append( num )
                else:
                    files[extn][stem]=[num,]
            else:
                files[extn]={stem:[num,]}
        else:
            names.append( f.name )
    # sort the numbers
    for extn in files:
        for stem in files[extn]:
            try:
                dsu = [ (int(s),s) for s in files[extn][stem] ]
            except:
                print(files[extn][stem] )
                raise
            dsu.sort()
            files[extn][stem] = [s[1] for s in dsu]
    return direcs, files, names

def print_a_folder( folder ):
    direcs , files, names = list_a_folder( folder )
    print("In folder:",folder)
    if len(direcs)>0:
        print("  Directories:")
        for d in sorted(direcs):
            print("\t",d)
    if len(names)>0:
        print("  Files:")
        for n in sorted(names):
            print("\t",n)
    if len(files)>0:
        for extn in files:
            for stem in files[extn]:
                nums = files[extn][stem]
                if len(nums)>3:
                    print("\t",stem+nums[0]+extn)
                    print("\t",stem+("?"*len(nums[0]))+extn,
                          " ... skip  ",len(nums)-2,"files")
                    print("\t",stem+nums[-1]+extn,
                          " ... total ",len(nums),"files")
                else:
                    for num in nums:
                        print("\t",stem+num+extn)
            

if __name__=="__main__":
    if len(sys.argv)>1:
        print_a_folder( sys.argv[1] )
    else:
        print_a_folder( "." )
                
                
        
                
