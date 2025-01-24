
import numpy as np

def norm( v ):
    modv = np.sqrt((np.array(v)*v).sum())
    return v / modv

def project( v ):
    """ project vector onto a pole figure """
    vm = norm( v )
    x = vm[0]/(vm[2] + 1)
    y = vm[1]/(vm[2] + 1)    
    return x,y

def ubi2pf( ubi, v = (0,0,1) ):
    """
    Get the RGB for an inverse pole figure
    """
    hkl = abs( norm( np.dot( ubi, v ) ))
    hkl.sort()
    return project( hkl )

REDPOS = project( (0,0,1 ) )
GREENPOS = project( (0,1,1) )
BLUEPOS = project( (1,1,1) )
CX = (REDPOS[0] + GREENPOS[0] + BLUEPOS[0])/3
CY = (REDPOS[1] + GREENPOS[1] + BLUEPOS[1])/3
dr0 = (REDPOS[0]   - CX)**2 + (REDPOS[1]   - CY)**2
dg0 = (GREENPOS[0] - CX)**2 + (GREENPOS[1] - CY)**2
db0 = (BLUEPOS[0]  - CX)**2 + (BLUEPOS[1]  - CY)**2


def pf2rgb( x, y ):   
    dr = np.arctan2(dr0, ((REDPOS[0]   - x)**2 + (REDPOS[1]   - y)**2))
    dg = np.arctan2(dg0, ((GREENPOS[0] - x)**2 + (GREENPOS[1] - y)**2))
    db = np.arctan2(db0, ((BLUEPOS[0]  - x)**2 + (BLUEPOS[1]  - y)**2))
    m = np.where( dr > dg, dr, dg )
    m = np.where( db > m, db, m )
    rgb = np.array( [dr, dg, db] )
    rgb = rgb / m
    return (rgb*255).astype(np.uint8)
            
def ubi2rgb( ubi ):
    x,y = ubi2pf( ubi )
    return pf2rgb( x, y )
