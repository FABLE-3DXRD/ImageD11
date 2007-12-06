

import unittest, os, sqlite3

from ImageD11 import colfile2db, columnfile

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
        colfile2db.colfile2db( "test.flt" , "test.db" )
        con = sqlite3.connect("test.db")
        cur = con.cursor()
        cur.execute("select * from peaks")
        i = 0
        for row in cur:
            for item in row:    
               i+=1
               self.assertEqual(i, item)
               
if __name__=="__main__":
    unittest.main()               
    
