
from __future__ import print_function


#
# Look up where this comes from.
# It was copied from somewhere and modified slightly


try:
    import Tkinter as Tk
except:
    import tkinter as Tk



class listdialog(Tk.Toplevel):
    """
    Dialog box for setting detector parameters
    Takes a list of strings and numbers
    """
    def __init__(self, parent, title = None, items=None, logic = None):
        Tk.Toplevel.__init__(self, parent)
        self.transient(parent)
        if title:
            self.title(title)
        self.logic = logic
        self.logicvars = {}
        self.parent = parent
        self.result = items
        body = Tk.Frame(self)
        self.initial_focus = self.body(body,items,logic)
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
        
    def body(self, master, items, logic=None):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        self.e=[]
        if items is not None:
            i=0
            keys=list(items.keys())
            keys.sort()
            self.keys=keys
            for key in keys:
                Tk.Label(master,text=key).grid(row=i)
                el=Tk.Entry(master)
                el.insert(Tk.END,items[key])
                el.grid(row=i,column=1)
                self.e.append(el)
                if logic is not None and key in logic:
                    val = logic[key]
                    self.logicvars[key] = Tk.IntVar()
                    self.logicvars[key].set(val)
                    b=Tk.Checkbutton(master,text="Vary?",
                                     variable=self.logicvars[key])
                    b.grid(row=i,column=2)
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
        self.fv = {}
        for item in self.e:
            k = self.keys[i]
            retdict[k]=item.get()
            if self.logic is not None and k in self.logic:
                self.fv[k]=self.logicvars[k].get()
            i=i+1
        self.result=retdict
        print(self.result)


class columnchooser(listdialog):
    """
    Dialog box for setting detector parameters
    Takes a list of strings and numbers
    """
    def __init__(self, parent, items, title="Choose two columns"):
        Tk.Toplevel.__init__(self, parent)
        self.transient(parent)
        if title:
            self.title(title)
        body = Tk.Frame(self)
        listbox1 = Tk.Listbox(body)
        listbox2 = Tk.Listbox(body)
        for i in items:
            listbox1.insert(Tk.END,i)
            listbox2.insert(Tk.END,i)
        body.pack(padx=5, pady=5)
        listbox1.pack()
        listbox2.pack()
        self.initial_focus=body
        self.buttonbox()
        self.grab_set()
        self.protocol("WM_DELETE_WINDOW", self.cancel)
        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))
        self.initial_focus.focus_set()
        self.wait_window(self)
        
