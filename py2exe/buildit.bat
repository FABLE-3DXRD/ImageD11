
d:
cd D:\wright\eclipse_workspaces\fabulous\fabio\trunk
c:\python24\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\xfab
c:\python24\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\ImageD11\trunk\
c:\python24\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\fabian
c:\python24\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\fabric
c:\python24\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\simul_farfiled_python
c:\python24\python setup.py build --compiler=mingw32 install

cd D:\wright\eclipse_workspaces\fabulous\ImageD11\trunk\py2exe
copy OpenGL.__init__.py c:\python24\lib\site-packages\OpenGl\__init__.py
copy OpenGL.Tk.__init__.py c:\python24\lib\site-packages\OpenGl\Tk\__init__.py

rmdir /s /q dist
c:\python24\python py2exesetup.py py2exe
copy C:\WINDOWS\system32\MSVCP71.dll dist\
copy c:\python24\lib\site-packages\OpenGl\Tk\win32-tk8.4\* dist\

