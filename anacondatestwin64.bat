
set SRC="\users\jon\Documents\fable\ImageD11"
set PYT=c:\anaconda\python

GOTO TEST


cd %SRC%
%PYT% setup.py build install 

%PYT% setup.py build sdist --formats=gztar,zip bdist_wininst

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


cd %SRC%\test\quantix
%PYT% testfitfilt.py

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

cd %SRC%/test
%PYT% test_put_incr.py

cd %SRC%\test\demo
%PYT% latred_new.py
%PYT% test.py

cd %SRC%/test/test_index_unknown
%PYT% test_index_unknown.py

cd %SRC%/test/test_mmo
%PYT% make_test_data.py

cd %SRC%/test
%PYT% testscale.py

cd %SRC%/test
%PYT% testcolumnfile.py

cd %SRC%/test/testcolfile2db
%PYT% testcolfile2db.py

cd %SRC%/fsrc/
%PYT% tst.py ..\test\nac_demo\peaks.out_merge_t200 ..\test\nac_demo\nac.prm
%PYT% test_xlylzl.py ..\test\nac_demo\peaks.out_merge_t200 ..\test\nac_demo\nac.prm

cd %SRC%/test/makemap
%PYT% ..\..\scripts\makemap.py -f test.flt -p test.prm -u map.ubi -t 0.02 -U map.map

%PYT% -c "import ImageD11; print ImageD11.__version__"

%PYT% -V

REM Icon added to ImageD11_gui.exe by using XN_resource editor
REM http://www.wilsonc.demon.co.uk/delphi.htm

cd %SRC%

:END

