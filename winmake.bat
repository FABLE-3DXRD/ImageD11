
set SRC=%~dp0
set PYT=python
cd %SRC%
%PYT% setup.py build --compiler=mingw32
%PYT% setup.py build bdist_wheel
cd dist
pip install ImageD11-1.6.0-cp27-cp27m-win_amd64.whl --no-deps -U

GOTO TEST




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


%PYT% -c "import ImageD11; print ImageD11.__version__"

%PYT% -V

cd %SRC%

:END

