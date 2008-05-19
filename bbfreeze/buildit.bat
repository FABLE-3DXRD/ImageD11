
REM Update all packages

d:
cd D:\wright\eclipse_workspaces\fabulous\fabio\trunk
svn update
c:\python25\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\xfab
svn update
c:\python25\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\ImageD11\trunk\
svn update
c:\python25\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\fabian
svn update
c:\python25\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\fabric
svn update
c:\python25\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\simul_farfiled_python
svn update
c:\python25\python setup.py build --compiler=mingw32 install
cd D:\wright\eclipse_workspaces\fabulous\PolyXSim
svn update
c:\python25\python setup.py build --compiler=mingw32 install

cd D:\wright\eclipse_workspaces\fabulous\ImageD11\trunk\bbfreeze
c:\python25\python freezem.py
