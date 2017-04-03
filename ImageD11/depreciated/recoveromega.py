#!/usr/bin/python


# ImageD11_v1.0 Software for beamline ID11
# Copyright (C) 2005-2007  Jon Wright
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

import sys

try:
    logfile = open(sys.argv[1],"r")
    pksfile = open(sys.argv[2],"r")
    newpksfile = open(sys.argv[3],"w")
except:
    print "Usage: %s logfile pksfile newpksfile"%(sys.argv[0])
    print " recovers omega angles from spec log files to fill information into"
    print " peaksearch output files"
    sys.exit()

# Read all of the lines in logfile into a dictionary of filenames/omega angles

lookups = {}

for line in logfile.readlines():
    if line.find(".edf")>0:
        try:
            # 12.34 /data/opid11/external/me001/mysample/mysample0012.edf
            # split()[0] = 12.34
            # split()[1] = /data/opid11/external/me001/mysample/mysample0012.edf
            # name = mysample0012.edf
            # Finally lookups['mysample0012.edf']="12.34" - voila.
            om = line.split()[0]
            fullname=line.split()[1]
            name=fullname.split('/')[-1]
            lookups[name]=om
        except:
            print line
            raise

logfile.close()

for line in pksfile.readlines():
    # Reading in the pksfile from the peaksearching script
    newpksfile.write(line) # echo line
    if line.find("File ")>0: # We found a filename
        name=line.split()[-1]
        try:
            newpksfile.write("# Omega = %s\n"%(lookups[name]) ) # So add the angle
        except:
            newpksfile.write("# Omega = -999\n" ) # So add the angle
# all done
pksfile.close()
newpksfile.close()
