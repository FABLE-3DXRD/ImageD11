
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


import Queue, threading
global stop_now
stop_now = False

class ImageD11_thread(threading.Thread):
    """ Add a stopping mechanism for unhandled exceptions """
    def __init__(self, myname="ImageD11_thread"):
        self.myname=myname
        threading.Thread.__init__(self)
    def run(self):
        global stop_now
        try:
            self.ImageD11_run()
        except:
            stop_now = True
            raise
    def ImageD11_stop_now(self):
        global stop_now
        if stop_now:
            print "Got a stop in",self.myname
        return stop_now
