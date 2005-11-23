#!/usr/bin/python

import sys

#!/bliss/users/blissadm/bin/python

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

def getheader(filename):
    """
    Reads a header from an edf file in 1024 byte chunks.
    Assumes enclosing { }
    Returns string enclosed
    Adds a filename key at the top
    """
    h="filename = "
    fp=open(filename,"rb")
    h=h+filename+";\n"
    s=fp.read(1024)
    if s.find("{")==-1:
        raise Exception("Not an edf file") 
    while 1:
        if s.find("}")>=0:
            h=h+s[0:s.find("}")+2]
            break
        else:
            h=h+s
        s=fp.read(1024)
    return h


if __name__=="__main__":
    keys=[]
    for arg in sys.argv[1:]:
        if arg.find("key=")>=0:
            keys.append(arg.split("=")[-1])
            continue
        h=getheader(arg)
        sys.stdout.write(arg+" ")
        if len(keys)>0:
            for key in keys:
                item=h[h.find(key):].split(";")[0]
                sys.stdout.write(" "+item+" ")
        else:
            sys.stdout.write(h)
        sys.stdout.write("\n")
                       
