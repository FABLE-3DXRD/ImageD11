


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

import sys

try:
   logfile = open(sys.argv[1],"r")
   pksfile = open(sys.argv[2],"r")
   newpksfile = open(sys.argv[3],"w")
except:
   print "Usage: %s logfile pksfile newpksfile"%(sys.argv[0])
   sys.exit()

# Read all of the scans into a dictionary of filenames/omega angles

lookups = {}

for line in logfile.readlines():
   if line.find(".edf")>0:
      try:
         om = line.split()[0]
         fullname=line.split()[1]
         name=fullname.split('/')[-1]
         lookups[name]=om
      except:
         print line
         raise

logfile.close()

for line in pksfile.readlines():
   newpksfile.write(line)
   if line.find("File ")>0:
      name=line.split()[-1]
      
      newpksfile.write("# Omega = %s\n"%(lookups[name]) )


pksfile.close()
newpksfile.close()
