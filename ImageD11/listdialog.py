

#
# Look up where this comes from.
# It was copied from somewhere and modified slightly



# from Numeric import *
import Tkinter as Tk



class listdialog(Tk.Toplevel):
    """
    Dialog box for setting detector parameters
    Takes a list of strings and numbers
    """
    def __init__(self, parent, title = None, items=None):
        Tk.Toplevel.__init__(self, parent)
        self.transient(parent)
        if title:
            self.title(title)
        self.parent = parent
        self.result = items
        body = Tk.Frame(self)
        self.initial_focus = self.body(body,items)
        body.pack(padx=5, pady=5)
        self.buttonbox()
        self.grab_set()
        if not self.initial_focus:
            self.initial_focus = self
        self.protocol("WM_DELETE_WINDOW", self.cancel)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))
        self.initial_focus.focus_set()
        self.wait_window(self)
    def body(self, master, items):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        self.e=[]
        if items!=None:
           i=0
           keys=items.keys()
           keys.sort()
           self.keys=keys
           for key in keys:
              Tk.Label(master,text=key).grid(row=i)
              el=Tk.Entry(master)
              el.insert(Tk.END,items[key])
              el.grid(row=i,column=1)
              self.e.append(el)
              i=i+1
           return self.e[0]
        
    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons
        box = Tk.Frame(self)
        w = Tk.Button(box, text="OK", width=10, command=self.ok, default=Tk.ACTIVE)
        w.pack(side=Tk.LEFT, padx=5, pady=5)
        w = Tk.Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=Tk.LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack()
    #
    # standard button semantics
    def ok(self, event=None):
        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return
        self.withdraw()
        self.update_idletasks()
        self.apply()
        self.cancel()
        
    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()
    #
    # command hooks
    def validate(self):
        return 1 # override
     
    def apply(self):
        retdict={}
        i=0
        for item in self.e:
           retdict[self.keys[i]]=item.get()
           i=i+1
        self.result=retdict
        print self.result

