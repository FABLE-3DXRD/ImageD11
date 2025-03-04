
from __future__ import print_function

"""
Row/column array type

Inherits from numpy.ndarray

To be used for lattice reduction to apply to g-vectors and patterson
peaks equally well and in a coherent way.
"""

from numpy import dot, array, allclose, asarray, fabs,\
    argmin, argmax, sqrt, argsort, take, sum, where, ndarray, eye,\
    zeros, cross
from numpy.linalg import inv, LinAlgError

# Confirm that dot'ting a 3x3 matrix with a 3x10 gives a 3x10
assert dot(eye(3), zeros( (3, 10) ) ).shape == (3, 10), \
    "Numpy dot insanity problem"
           
# It is unclear why it is not a 10x3 result (row/col vectors)

try:
    dot(eye(3), zeros( (10, 3) ) )
    raise Exception("Numpy dot insanity problem")
except ValueError: 
    pass
except:
    print("Unexpected exception when checking numpy behaviour")
    raise



# Based on http://www.scipy.org/Subclasses
class rc_array(ndarray):
    """ 
    Row/Column array
    Represent a list of row or column vectors
    """
    def __new__(subtype, data, direction=None, dtype=None, copy=False):
        """ 
        Mostly as from example
        direction is one of row / column
        """
        subarr = array(data, dtype=dtype, copy=copy)
        subarr = subarr.view(subtype)
        if direction is not None:
            subarr.direction = direction
        elif hasattr(data, 'info'):
            subarr.direction = data.direction
        return subarr
    
    def __array_finalize__(self, obj):
        """ 
        Fill in a default row arg to direction
        self/obj??
        """
        self.direction = getattr(obj, 'direction', 'row' )
        

    def __str__(self):
        """ Used for printing """
        desc = \
"""{%(data)s in %(direction)s direction}"""
        return desc % { 'data': super.__str__(self), 
                        'direction' : self.direction }

    def __iter__(self):
        """
        Iterate depending on rows columns
        Use to get [ v for v in rr_array ]
        """
        # print "iter called"
        if self.direction == 'row':
            return ndarray.__iter__(self)
        elif self.direction == 'col':
            return ndarray.__iter__(self.T)
        else:
            raise Exception("rc_array with direction not in row|col")

    def norm2(self):
        """ sum(v*v,axis=? for row or col """
        return sum( self*self, axis=self.vector_axis())

    def vector_axis(self):
        """ The axis which has the 3 on it """
        return ['col', 'row'].index(self.direction)

    def nb_vector_axis(self):
        """ The axis which has the n on it """
        return ['row', 'col'].index(self.direction)

    def other_direction(self):
        """ The one which is not self.direction
        """
        return ['row', 'col'][self.vector_axis()]

    def check(self):
        """ 
        Ensure we have an rc_array which is well behaved 
        Pattern assert(check(v)) should disappear in optimiser
        """
        assert hasattr(self, 'direction')
        assert self.direction in ['row', 'col']
        if len(self.shape) == 1:
            assert self.shape[0] == 3, str(self)
        elif len(self.shape) == 2:
            if self.direction == 'row':
                assert self.shape[1] == 3
            if self.direction == 'col':
                assert self.shape[0] == 3
        else:
            raise Exception("Only 1D or 2D rc_arrays allowed so far")
        return True

    def nvectors(self):
        return self.shape[self.nb_vector_axis()]

    def flip(self, mat):
        """
        Flip row to col or col.row
        """
        assert self.check()
        assert mat.shape == (3, 3)
        if self.direction == 'row':
            ret = dot( mat , self.T)
        if self.direction == 'col':
            ret = dot( mat, self).T
        ret = rc_array( ret, direction = self.other_direction())
        ret.check()
        assert ret.shape == self.shape[::-1],"Shape mismatch in flip"
        return ret

    def inv(self):
        """
        Inverse matrix of self
        """
        assert self.shape == (3,3)
        ret = inv(self)
        return rc_array(ret, self.other_direction())

if __name__=="__main__":


    v = rc_array([1,2,3],direction='row')
    assert v.other_direction() == 'col'

    v = rc_array([[1,2,3],[4,5,6]],direction='row')

    assert v.other_direction() == 'col'
    assert v.flip(eye(3)).direction == 'col'
    assert allclose( v.flip(eye(3)) , v.T )
    assert v.norm2().shape == ( v.shape[ v.nb_vector_axis() ] ,)
    assert allclose( v.norm2(), array([ 1+4+9, 16+25+36 ])), str(v)
    

    v = rc_array(array( [[1,2,3],[4,5,6]] ).T ,direction='col')

    assert v.other_direction() == 'row'
    assert v.flip(eye(3)).direction == 'row'
    assert allclose( v.flip(eye(3)) , v.T )

    assert v.norm2().shape == ( v.shape[ v.nb_vector_axis() ] ,)

    assert allclose( v.norm2(), array([ 1+4+9, 16+25+36 ])), str(v)

