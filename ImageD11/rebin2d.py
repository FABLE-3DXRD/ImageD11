


# ImageD11_v0.4 Software for beamline ID11
# Copyright (C) 2005  Jon Wright
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


"""
Function for rebinning in two dimensions

Expected to be very slow and then re-implemented in C once it works, correctly

"""

from math import ceil,floor,sqrt


class polygon:
    """
    Represents a 2D polygon
    """
    def __init__(self,xy):
        """
        xy are a list of pairs of xy points (directed)
        """
        self.xy=xy
    def walkaroundintegervertices(self):
        """
        Generate a list of points along the edges of integers
        """
        path=[]
        l=[ item for item in self.xy ]
        l.append(self.xy[0])
        for j in range(len(l)-1):
            p1 = l[j]
            p2 = l[j+1]
            path.append(l[j])
            # Find points along edge
            intersects=[]
            for i in range(ceil(p1[0]),ceil(p2[0])):
                # Intersections on zeroth axis
                g = 1.*(p2[1]-p1[1])/(p2[0]-p1[0])
                point=[i,p1[1]+g*(i-p1[0])]
                intersects.append( [distance(p1,point),point] )
            for i in range(ceil(p1[1]),ceil(p2[1])):
                # Intersections on oneth axis
                g = 1.*(p2[0]-p1[0])/(p2[1]-p1[1])
                point=[p1[0]+g*(i-p1[1]),i]
                intersects.append( [distance(p1,point),point] )
            for i in range(ceil(p2[0]),ceil(p1[0])):
                # Intersections on zeroth axis
                g = 1.*(p2[1]-p1[1])/(p2[0]-p1[0])
                point=[i,p1[1]+g*(i-p1[0])]
                intersects.append( [distance(p1,point),point] )
            for i in range(ceil(p2[1]),ceil(p1[1])):
                # Intersections on oneth axis
                g = 1.*(p2[0]-p1[0])/(p2[1]-p1[1])
                point=[p1[0]+g*(i-p1[1]),i]
                intersects.append( [distance(p1,point),point] )
            if len(intersects)>0:
                intersects.sort()
                # print "intersects",intersects
                for d,point in intersects:
                    path.append(point)

        self.path=path
                
                    

def distance(p1,p2):
    return sqrt(p1[0]*p1[0]+p2[0]*p2[0])


def main():
    print "started"
    vertices = [ [ 1.1, 0.9] ,
                 [ 3.2, 1.1] ,
                 [ 3.2, 4.2] ,
                 [ 1.2, 3.1] ]
    obj = polygon(vertices)
    obj.walkaroundintegervertices()

if __name__=="__main__":
    main()

