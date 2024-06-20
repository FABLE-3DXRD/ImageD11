import unittest

""" Write some test cases for the columnfile stuff """

from ImageD11.columnfile import colfile_to_hdf, PandasColumnfile as columnfile

class testgeom(unittest.TestCase):
    def setUp(self):
        with open("testgeom.flt", "w") as f:
            f.write("""
# chi = 0.0
# distance = 5881.14814327
# o11 = 1
# o12 = 0
# o21 = 0
# o22 = 1
# omegasign = 1.0
# t_x = 31.8613635188
# t_y = -28.5519205159
# t_z = 0.0
# tilt_x = -0.00319993715326
# tilt_y = -0.00227706383676
# tilt_z = 0.0172853900689
# wavelength = 0.2469981
# wedge = 0.0761143411257
# y_center = 1080.16471254
# y_size = 1.5
# z_center = 978.22872971
# z_size = 1.5
#  sc  fc  omega  Number_of_pixels
20  30  0  1
21  31  30  1
22  32  60  1
23  33  70  1
24  34  80  1
""")

    def test1(self):
        c = columnfile("testgeom.flt")
        c.updateGeometry()

    def testfilter(self):
        c = columnfile("testgeom.flt")
        d = c.copy()
        d.filter(abs(d.sc - 22) > .1)
        self.assertTrue(d.nrows == 4)
        d.filter(abs(d.sc - 22) > .1)
        self.assertTrue(d.nrows == 4)
        self.assertTrue(c.nrows == 5)
        d = c.copy()
        self.assertTrue(d.nrows == 5)

    def testfilter2(self):
        # What's with the bigarray call?
        c = columnfile("testgeom.flt")
        d = c.copy()
        # a = d.bigarray
        d.filter(abs(d.sc - 22) > .1)
        self.assertTrue(d.nrows == 4)
        d.filter(abs(d.sc - 22) > .1)
        self.assertTrue(d.nrows == 4)
        self.assertTrue(c.nrows == 5)
        d = c.copy()
        self.assertTrue(d.nrows == 5)


class test1(unittest.TestCase):
    def setUp(self):
        f = open("test.col", "w")
        f.write("""#  x  a  b  c
0 1 2 3
5 6 1 2
0 4 2 3
32 5 6 2
33 7 9 0
""")
        f.close()

    def test_to_hdf(self):
        c = columnfile("test.col")
        colfile_to_hdf(c, "testcolfile.hdf")
        h = columnfile("testcolfile.hdf")
        for t in c.titles:
            assert ((c.getcolumn(t) == h.getcolumn(t)).all())
            assert t in h.titles

    def testaddcol1(self):
        c = columnfile("test.col")
        c.addcolumn([5, 4, 1, 2, 0], 'alice')
        self.assertEqual(c.titles[-1], 'alice')
        self.assertEqual(list(c.alice), [5, 4, 1, 2, 0])
        self.assertEqual(c.ncols, 5)
        self.assertEqual(c.nrows, 5)

    def testaddcol2(self):
        c = columnfile("test.col")
        c.addcolumn([5, 4, 1, 2, 9], 'a')
        self.assertEqual(c.titles.index("a"), 1)
        self.assertEqual(list(c.a), [5, 4, 1, 2, 9])
        self.assertEqual(c.ncols, 4)
        self.assertEqual(c.nrows, 5)

    def test1(self):
        c = columnfile("test.col")
        self.assertEqual(list(c.b), [2, 1, 2, 6, 9])
        c.removerows("b", [1])
        self.assertEqual(c.nrows, 4)
        self.assertEqual(list(c.b), [2, 2, 6, 9])

    def test2(self):
        c = columnfile("test.col")
        assert list(c.b) == [2, 1, 2, 6, 9]
        c.removerows("b", [2])
        assert c.nrows == 3
        assert list(c.b) == [1, 6, 9]

    def test3(self):
        c = columnfile("test.col")
        assert list(c.a) == [1, 6, 4, 5, 7]
        assert list(c.b) == [2, 1, 2, 6, 9]
        c.removerows("a", [1.001], tol=0.1)
        c.removerows("a", [4])
        assert c.nrows == 3
        assert list(c.b) == [1, 6, 9]

    def test4(self):
        c = columnfile("test.col")
        assert list(c.a) == [1, 6, 4, 5, 7]
        c.removerows("a", [1, 4])
        assert c.nrows == 3
        assert list(c.b) == [1, 6, 9]

    def testsort(self):
        c = columnfile("test.col")
        c.sortby('a')
        assert (list(c.a.astype(int)) == [1, 4, 5, 6, 7])

    def testreorder(self):
        c = columnfile("test.col")
        c.reorder([4, 3, 2, 1, 0])
        assert (list(c.a) == [7, 5, 4, 6, 1])

    def testmissingfile(self):
        e = False
        try:
            columnfile("a_filename_that_does_not_exist.flt")
        except:
            # got an exception
            e = True
        assert (e)


if __name__ == '__main__':
    unittest.main()
