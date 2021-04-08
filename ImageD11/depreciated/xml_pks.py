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
# MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA



"""
Script for peaksearching images from the command line

Uses the connectedpixels extension for finding blobs above a threshold
and the blobcorrector(+splines) for correcting them for spatial distortion

Defines one function (peaksearch) which might be reused
"""

import xml.etree.ElementTree as ET

class p:
    root = ET.Element("root")
    def add_option(self, *args, **kwds):
        junk = []
        if len(args)==2:
            short, int = args
            e = ET.Element("option")
            e.text = int[2:] # strip --
            junk.append( ET.SubElement(e, "short_name") )
            junk[-1].text = short[1:] # strip -
        elif len(args) == 1:
            long = args[0]
            e = ET.Element("option")
            e.text = int[2:] # strip --
        else:
            raise Exception("Not one or two args")
        keys = list(kwds.keys())
        for k in keys:
            junk.append( ET.SubElement(e, k) )
            junk[-1].text = str(kwds[k])
        self.root.append(e)

    def save(self, filename):
        self.indent(self.root, 1)
        tree = ET.ElementTree(self.root)
        tree.write(filename)


    def indent(self, elem, level=0):
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            for elem in elem:
                self.indent(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i


if __name__=="__main__":
        parser = p()
        import ImageD11.peaksearcher
        newparser = ImageD11.peaksearcher.get_options(parser)
        newparser.save("peaksearch.xml")

