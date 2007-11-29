
set SRC="\wright\eclipse_workspaces\fabulous\imaged11\trunk"
set PYT=c:\python24\python
cd %SRC%
%PYT% setup.py build --compiler=mingw32 install 

%PYT% setup.py build --compiler=mingw32 sdist --formats=gztar,zip bdist_wininst
GOTO TEST

REM Py2exe is not quite automatic. Doubtless this is version and computer specific...
cd py2exe
copy OpenGL.__init__.py c:\python24\lib\site-packages\OpenGl\__init__.py
copy OpenGL.Tk.__init__.py c:\python24\lib\site-packages\OpenGl\Tk\__init__.py
del /q dist build
%PYT% py2exesetup.py py2exe
copy c:\python24\lib\site-packages\OpenGl\Tk\win32-tk8.4\* dist\
cd dist
REM del /q _na* _ns* ssl.pyd

:TEST

cd %SRC%\test\demo
%PYT% test.py

cd %SRC%\test\quantix
%PYT% fitgrain.py g3.pars g3.ubi g3.flt new.pars

cd %SRC%\test
%PYT% test_peaksearch.py

cd %SRC%\test
## takes really ages
REM python test_peaksearch.py ALL

cd %SRC%\test\gv_general
%PYT% test_gv_general.py

cd %SRC%\test\testconnectedpixels
%PYT% testconnectedpixels.py

cd %SRC%\test\testlabelimage
%PYT% testlabelimage.py


%PYT% -c "import ImageD11; print ImageD11.__version__"

%PYT% -V

REM Icon added to ImageD11_gui.exe by using XN_resource editor
REM http://www.wilsonc.demon.co.uk/delphi.htm

cd %SRC%

:END

