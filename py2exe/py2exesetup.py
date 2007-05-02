


from distutils.core import setup
import py2exe

from distutils.filelist import findall
import os
import matplotlib
matplotlibdatadir = matplotlib.get_data_path()
matplotlibdata = findall(matplotlibdatadir)
matplotlibdata_files = []
for f in matplotlibdata:
    dirname = os.path.join('matplotlibdata', f[len(matplotlibdatadir)+1:])
    matplotlibdata_files.append((os.path.split(dirname)[0], [f]))



        
setup(
     console = ["../scripts/ImageD11_gui.py","../scripts/peaksearch.py",
                "../scripts/ImageD11Server.py" ],
     options={
             'py2exe': {
                        'excludes' : ['numpy','wx'],
                        'includes' : ['OpenGL.GL', 'OpenGL.Tk'],
                        'packages' : ['matplotlib', 'pytz'],
                        'dll_excludes':['libgdk-win32-2.0-0.dll',
                             'libgobject-2.0-0.dll',
                             'libgdk_pixbuf-2.0-0.dll',
                             'wxmsw26uh_vc.dll']
                       }
            },
     data_files=matplotlibdata_files
)

# Sorry - no idea where the dll excludes come from in the first place
# ... and copy pyopengl directory to your dist directory
print
print "Have you patched your pyopengl and copied togl.dll to dist "
print "eg:"
print "copy OpenGL.__init__.py c:\\python24\\lib\\site-packages\\OpenGl\\__init__.py"
print "copy OpenGL.Tk__init__.py c:\\python24\\lib\\site-packages\\OpenGl\\Tk\\__init__.py"
print "copy c:\\python24\\lib\\site-packages\\OpenGl\\Tk\\win32-tk8.4\\* dist\\"

