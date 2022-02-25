from __future__ import print_function

import unittest, os, sqlite3

from ImageD11 import columnfile

class t1(unittest.TestCase):
    def setUp(self):
        open("test.flt","w").write(
"""# col1 col2 col3
1 2 3
4 5 6
7 8 9
""")
        try:
            os.remove("test.db")
        except:
            pass

    def test1(self):
        """ vague check
        """
        columnfile.colfile2db( "test.flt" , "test.db" )
        con = sqlite3.connect("test.db")
        cur = con.cursor()
        cur.execute("select * from peaks")
        i = 0
        for row in cur:
            for item in row:    
               i+=1
               self.assertEqual(i, item)
import os, time

class t2(unittest.TestCase):
    def setUp(self):
        if os.path.exists("nac.db"):
            return
        start = time.time()
        columnfile.colfile2db( os.path.join( os.path.split(__file__)[0],
            "..","nac_demo","peaks.out_merge_t200") , "nac.db" )
        print( "write db",time.time()-start)

    def test1(self):
        # read back
        start = time.time()
        con = sqlite3.connect("nac.db")
        cur = con.cursor()
        cur.execute("select * from peaks")
        dat = [ [ v for v in row] for row in cur]
        print( "read db",time.time()-start, len(dat[0]), len(dat))
        

if __name__=="__main__":
    unittest.main()               
    
