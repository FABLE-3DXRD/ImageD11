
# coding: utf-8

# In[1]:

from ImageD11 import unitcell, indexing
import numpy as np
import matplotlib.pylab as plt
#get_ipython().magic(u'matplotlib inline')


# In[2]:

u = unitcell.unitcell([3.,3.,3.1,90,90,90],'P')
u.B[1,0]=u.B[2,1]=u.B[0,2]=0
u.B[0,1]=u.B[1,2]=u.B[2,0]=0


# In[3]:

u.makerings(3.11/3)


# In[4]:

rng = range(-2,3)
hkls = np.array([(h,k,l) for h in rng for k in rng for l in rng]).T


# In[5]:

UB0 = u.B
print UB0


# In[6]:

gve = np.dot( UB0, hkls )


# In[7]:

modg= np.sqrt((gve*gve).sum(axis=0))


# In[8]:

plt.plot(modg, np.arctan2(gve[1],gve[2]),"+")


# In[9]:

order = np.argsort(modg)
for i in order:
    print i,hkls.T[i],modg[i],gve.T[i]


# In[9]:




# In[10]:

u.makerings(0.6)
for d in u.ringds:
    print d,u.ringhkls[d][0]


# In[11]:

u.orient( 1, gve.T[87], 2,gve.T[68])


# In[12]:

len(u.UBIlist)


# In[13]:

# does it index the peaks gve.T[4] and gve.T[9]?
for ubi in u.UBIlist:
    print np.dot(ubi,gve.T[87]),np.dot(ubi,gve.T[68]),np.linalg.det(ubi)


# In[14]:

# how many peaks are indexed?
def countpks(ubi, gve, tol = 1e-6):
    hkls = np.dot(ubi,gve)
    ihkls = np.round(hkls)
    dh = hkls - ihkls
    drlv = np.sqrt((dh*dh).sum(axis=0))
    return (drlv<tol).sum()


# In[15]:

for ubi in u.UBIlist:
    print countpks(ubi,gve),"\n",ubi


# In[16]:

print u.UBIlist[0]

