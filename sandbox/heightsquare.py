
import numpy as np

def heightsquare(r,theta):
    """ Height of a projected square versus x """
    p = ((-1,-1,1,1),
         (1,-1,-1,1))
    s = np.sin(np.radians(theta))
    c = np.cos(radians(theta))
    m = ((c, -s),(s, c))
    rp = np.dot(m,p)/2
    order = np.argsort(rp[0])
    o = []
    xp = rp[0][order]
    yp = rp[1][order]
    for v in r:
        if v < xp[0]:
            o.append(0)
            continue
        if v > xp[-1]:
            o.append(0)
            continue
        if v < xp[1]:
            h1 = abs(yp[1]-yp[0])*(v-xp[0])/(xp[1]-xp[0])
            h2 = abs(yp[2]-yp[0])*(v-xp[0])/(xp[2]-xp[0])
            o.append( h1+h2 )
            continue
        if v > xp[2]:
            h1 = abs(yp[3]-yp[2])*(xp[3]-v)/(xp[3]-xp[2])
            h2 = abs(yp[3]-yp[1])*(xp[3]-v)/(xp[3]-xp[1])
            o.append( h1+h2 )
            continue
        h1 = abs(yp[1]-yp[0])
        h2 = abs(yp[2]-yp[0])*(xp[1]-xp[0])/(xp[2]-xp[0])
        o.append( h1 + h2 )
    return o

if __name__=="__main__":
    x=np.linspace(-1,1,2000)
    from pylab import *
    h = heightsquare
    plot(x,h(x,0),"-")
    plot(x,h(x,1),"-")
    plot(x,h(x,10),"-")
    plot(x,h(x,50),"-")
    plot(x,h(x,45),"-")
    plot(x,h(x,79),"-")
    plot(x,h(x,99),"-")
    show()
