

from ImageD11.sym_u import cubic
from ImageD11.columnfile import columnfile
import numpy as np, sys, xfab.tools

# Cubic symmetry operators
m = np.array(cubic().group)
# print m.shape
def process_file(f, g, ip, jp, npoint):
    """
    f = gff file / sinfit file
    g = columnfile with allgrains output
    i, j = origin position in sample
    n = point number
    """
    try:
        c = columnfile( f+".gff")
    except:
        # print "No grains"
        # print
        # print
        return
    u = []
    u_uniq = []
    for i in range(c.nrows):
        ut =  np.array( [[ c.U11[i], c.U12[i], c.U13[i] ],
                         [ c.U21[i], c.U22[i], c.U23[i] ],
                         [ c.U31[i], c.U32[i], c.U33[i] ]])
        uniq = np.dot(m, ut).trace(axis1=1,axis2=2).argmax()
        u_uniq.append( np.dot( m[uniq], ut ))
        u.append(ut)
    pos = False
    axes = []
    xs = []
    ys = []
    allints = []
    ints = []
    for line in open(f+".sinfit"):
        if line.find("Constant")>0:
            pos = True
            if len(ints)>0:
                x = np.array(ints)
                allints.append([ x.mean(),x.std(),x.max(),x.min(),
                                 np.median(x),len(x) ] )
            ints = []
            continue
        if line[0] == "#":
            continue
        if pos:
            pos = False
            axis, x, y = line.split()
            axes.append(axis)
            xs.append(x)
            ys.append(y)
            continue
        vals = line.split()
        if len(vals) == 4:
            ints.append(float(vals[-1]))
        
    if len(ints)>0:
        x = np.array(ints)
        allints.append([ x.mean(),x.std(),x.max(),x.min(),
                         np.median(x),len(x) ] )
    # print len(allints)
    for i in range(c.nrows):
        g.write("%d  %d  %.2f  %.2f  "%(npoint,i,ip,jp))
        # g.write( ("%.8f  "*9)%tuple(u[i].ravel()))
        g.write( ("%.6f  "*3)%tuple(xfab.tools.u_to_rod(u_uniq[i])))
        g.write( ("%.6f  "*9)%tuple(u_uniq[i].ravel()))
        g.write( ("%s  "*3)%(axes[i], xs[i] , ys[i]))
        g.write(("%.6g  "*5 + "%d  \n")%tuple(allints[i]) )



if __name__=="__main__":

    titles="#  npoint  ngrain  i  j  R1  R2  R3  U11  U12  U13  U21  U22  U23  " +\
           "U13  U23  U33  axis  x  y  Imean  Istd  Imax  Imin  Imed  Npks\n"

    if len(sys.argv)==2:
        sys.stdout.write(titles)
        process_file( sys.argv[1], sys.stdout, 0, 0, int(sys.argv[1]) )
        sys.exit()


    n = 0
    i=-12
    j=-12
    # g = open("allgrains.flt","w")
    g = sys.stdout
    g.write(titles)
    while i < 12.1:
        while j < 12.1:
            # print "# Point", i, j
            process_file( str(n) , g, i, j, n)
            n += 1
            j += 0.5
        i+= 0.5
        j = -12
        
