set LDFLAGS="-static-libgfortran -static-libgcc -static -lgomp -shared -O3"
f2py -m fImageD11 -c fImageD11.f90 --f90flags="-fopenmp" -lgomp -lpthread --fcompiler=gnu95 --compiler=mingw32

echo "next"
cd ..
python setup.py build
copy fsrc\fImageD11.pyd build\lib.win-amd64-2.7
set PYTHONPATH=%~dp0build\lib.win-amd64-2.7
cd fsrc
python tst.py ..\test\nac_demo\peaks.out_merge_t200 ..\test\nac_demo\nac.prm
