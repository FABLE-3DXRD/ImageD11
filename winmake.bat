d:
set SRC="\wright\eclipse_workspaces\fabulous\imaged11\trunk"
cd %SRC%
setup.py build --compiler=mingw32 install 

setup.py build --compiler=mingw32 sdist --formats=gztar,zip bdist_wininst
REM Py2exe is not quite automatic. Doubtless this is version and computer specific...
cd py2exe
copy OpenGL.__init__.py c:\python24\lib\site-packages\OpenGl\__init__.py
copy OpenGL.Tk.__init__.py c:\python24\lib\site-packages\OpenGl\Tk\__init__.py
del /q dist build
py2exesetup.py py2exe
copy c:\python24\lib\site-packages\OpenGl\Tk\win32-tk8.4\* dist\
cd dist
del /q _na* _ns* ssl.pyd


cd %SRC%\test\demo
python test.py

cd %SRC%\test\quantix
python fitgrain.py g3.pars g3.ubi g3.flt new.pars

cd %SRC%\test
python test_peaksearch.py

cd %SRC%\test
## takes really ages
REM python test_peaksearch.py ALL

cd %SRC%\test\gv_general
python test_gv_general.py

cd %SRC%\test\testconnectedpixels
python testconnectedpixels.py

cd %SRC%\test\testlabelimage
python testlabelimage.py


python -c "import ImageD11; print ImageD11.__version__"


REM Icon added to ImageD11_gui.exe by using XN_resource editor
REM http://www.wilsonc.demon.co.uk/delphi.htm

:END

