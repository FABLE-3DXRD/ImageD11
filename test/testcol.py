
from ImageD11 import columnfile

open("test.col","w").write("""#  x  a  b  c
0 1 2 3
5 6 1 2
0 4 2 3
32 5 6 2
""")

c  = columnfile.columnfile( "test.col" ) 
assert list(c.b) == [ 2, 1, 2, 6]
c.removerows( "b" , [1] )
assert c.nrows == 3
assert list(c.b) == [ 2, 2, 6]

c  = columnfile.columnfile( "test.col" ) 
assert list(c.b) == [ 2, 1, 2, 6]
c.removerows( "b" , [2] )
assert c.nrows == 2
assert list(c.b) == [ 1, 6]

c  = columnfile.columnfile( "test.col" ) 
assert list(c.a) == [ 1, 6, 4, 5]
c.removerows( "a" , [1] )
c.removerows( "a" , [4] )
assert c.nrows == 2
assert list(c.b) == [ 1, 6]

c  = columnfile.columnfile( "test.col" ) 
assert list(c.a) == [ 1, 6, 4, 5]
c.removerows( "a" , [1, 4] )
assert c.nrows == 2
assert list(c.b) == [ 1, 6]

print "Passed"
