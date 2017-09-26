

set SRC=%~dp0
set PYT=python
cd %SRC%
%PYT% setup.py build --compiler=mingw32
%PYT% setup.py build bdist_wheel
cd dist
pip install ImageD11-1.8.0-cp27-cp27m-win_amd64.whl --no-deps -U


GOTO TEST



:TEST

cd %SRC%\test\demo
%PYT% latred_new.py
%PYT% test.py

cd %SRC%\test
%PYT% run_tests.py


cd %SRC%\test\quantix
%PYT% testfitfilt.py

cd %SRC%\test
%PYT% test_peaksearch.py

cd %SRC%\test
## takes really ages
REM python test_peaksearch.py ALL

cd %SRC%\test\ken_simul
%PYT% idx.py


%PYT% -c "import ImageD11; print ImageD11.__version__"

%PYT% -V

cd %SRC%

:END

