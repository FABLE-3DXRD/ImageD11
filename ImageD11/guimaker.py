

# This is copied from a book somewhere
# either programming python or the cookbook.
# look it up and give credit please..! Probably the book from Mark Lutz??

import sys
from Tkinter import *

class GuiMaker(Frame):   # Inherit from Tk frame
   menuBar=[]
   def __init__(self,parent=None):
      Frame.__init__(self,parent)
      self.pack(expand=YES, fill=BOTH)
      self.statusbar=Label(self,text="Welcome to ImageD11_v0.4")
      self.statusbar.pack(side=BOTTOM)
      self.start()
      self.makeWidgets()
      self.makeMenuBar()

   def makeMenuBar(self):
      menubar = Menu(self.master)
      self.master.config(menu=menubar)
      for (name, key, items) in self.menuBar:
         pulldown = Menu(menubar)
         self.addMenuItems(pulldown,items)
         menubar.add_cascade(label=name, underline=key, menu=pulldown)

   def addMenuItems(self, menu, items):
      for item in items:
         if item=="separator":
            menu.add_separator({})
         else:
            menu.add_command(label = item[0], underline = item[1], command=item[2] )


