
"""
Script to put the menus together and build an appli with each bit
being relatively clutterfree

Each of guiindexer, guipeaksearch and guitransformer offer members
menuitems to put in menus. The overall gui (imaged11_gui) inherits
from this I think

# The code was copied from a book somewhere
# probably programming python or maybe cookbook.
# Most probably programming python by Mark Lutz
"""
try:
    import Tkinter as Tk
except:
    # python 3 ? 
    import tkinter as Tk

class GuiMaker(Tk.Frame):   # Inherit from Tk frame
    """
    You must inherit from this class and implement the start and makeWidgets
    methods
    """
    menuBar=[]
    def __init__(self,parent=None):
        Tk.Frame.__init__(self,parent)
        self.pack(expand=Tk.YES, fill=Tk.BOTH)
        import ImageD11
        self.statusbar=Tk.Label(self,text="Welcome to ImageD11 version "+ImageD11.__version__)
        self.statusbar.pack(side=Tk.BOTTOM)
        self.start()
        self.makeWidgets()
        self.makeMenuBar()

    def makeMenuBar(self):
        menubar = Tk.Menu(self.master,tearoff=0)
        self.master.config(menu=menubar)
        for (name, key, items) in self.menuBar:
            pulldown = Tk.Menu(menubar)
            self.addMenuItems(pulldown,items)
            menubar.add_cascade(label=name, underline=key, menu=pulldown)

    def addMenuItems(self, menu, items):
        for item in items:
            if item=="separator":
                menu.add_separator({})
            else:
                menu.add_command(label = item[0], underline = item[1], command=item[2] )
