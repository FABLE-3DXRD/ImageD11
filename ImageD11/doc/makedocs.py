

import pydoc, os, ImageD11,inspect, shutil
pydoc.writedoc("ImageD11")
shutil.copyfile("ImageD11.html","index.html")
done=[None]
for file in os.listdir(ImageD11.__path__[0]):
    path = os.path.join(ImageD11.__path__[0], file)
    modname = inspect.getmodulename(file)
    if modname not in done: # and modname != '__init__':
        pydoc.writedoc("ImageD11."+modname)
        done.append(modname)
