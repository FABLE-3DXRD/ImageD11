

"""
Copy an ImageD11 filtered peaks file into an sqlite3 database
"""

import sqlite3 as database_module
# Perhaps other modules follow the same api.
# Doubtless one does not get away with using a filename?

from ImageD11 import columnfile

def colfile2db( colfilename, dbname ):
    """
    Read the columnfile into a database
    Ignores parameter metadata (not used yet)
    """
    colf = columnfile.columnfile( colfilename )
    dbo = database_module.connect( dbname )
    c = dbo.cursor()
    # Build up columnames and types to make table
    tablecols = []
    # Not allowed for sql to have ^ in string
    colf.titles = [t.replace("^","_pow_") for t in colf.titles]
    for name in colf.titles:
        if name in columnfile.INTS:
            tablecols.append(name + " INTEGER")
            continue
        if name in columnfile.FLOATS:
            tablecols.append(name + " REAL")
            continue
        tablecols.append(name + " REAL")
    c.execute("create table peaks \n( " + \
              " , ".join(tablecols)     + " ) ; \n" )
    # Make a format string for inserting data
    ins = "insert into peaks values ("  + \
          ",".join(["?"]*colf.ncols)    +") ;"
    # insert the data
    for i in range(colf.nrows):
        c.execute( ins , tuple(colf.bigarray[:,i]) )
    c.close()
    dbo.commit()
    dbo.close()

if __name__=="__main__":
    import sys
    try:
        colfile2db(sys.argv[1], sys.argv[2])
    except:
        print "Usage: %s fltfile sqlite3_db_file"
        raise
                                       
    
