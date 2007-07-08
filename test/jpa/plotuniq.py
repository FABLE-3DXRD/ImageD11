from matplotlib.pylab import *
from matplotlib.colors import cnames
keys = cnames.keys()




def HSVtoRGB(  h, s, v ):
    """ http://www.cs.rit.edu/~ncs/color/t_convert.html """
    if s == 0 :
    	return v,v,v
    h /= 60;		#	 sector 0 to 5
    i = floor( h );
    f = h - i;		#	 factorial part of h
    p = v * ( 1 - s );
    q = v * ( 1 - s * f );
    t = v * ( 1 - s * ( 1 - f ) );

    if i == 0: return v,t,p	
    if i == 1: return q,v,p
    if i == 2: return p,v,t		
    if i == 3: return p,q,v
    if i ==4 : return t,p,v
    return v,p,q	








z = []
i = []
g=0
filloff = 0.
for line in open(sys.argv[1],"r").readlines():
    if line.find("ubi[")>0 and len(z)>0:
        a = argsort
        g+=12.
        c = HSVtoRGB(g%360,0.8,0.9)
        # (g*0.4)%1, (g*0.45)%1, (g*0.55)%1
        if max(c)>1.:
            print c,g
        if z[0]>filloff and max(i)>1.5e9:        
            cf = HSVtoRGB(g%360,0.3,0.95)

            fill([0]+i+[0],array([z[0]]+z+[z[-1]])-3.15,label="%d"%(g),facecolor=cf,
                    edgecolor=c)
            filloff = z[-1]
        plot(i,array(z)-3.15,marker="+",linestyle="-", color = c)
        z=[]
        i=[]
    if line.find("ubi[")>0:
        z=[]
        i=[]
    if line.find("520\\scan")==-1:
        continue
    items = line.split()
    if int(items[0])<15:
        continue
    zh = float(items[1])
    if len(z) > 0 and zh > z[-1]:
        z.append(float(items[1]))
        diff = zh - z[0]
        i.append(float(items[-1]))
    else:
        if len(z) == 0:
            z.append(float(items[1]))
            diff = zh - z[0]
            i.append(float(items[-1]))
        else:
            print line.rstrip(),z,zh,len(z)
#legend(loc=0)
title("Diffracted intensity versus height")
xlabel("Total intensity")
ylabel("Height in sample")
show()
