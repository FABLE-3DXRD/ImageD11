
cd "\documents and settings\wright\eclipse_new\imaged11"
setup.py build --compiler=mingw32 install 
cd test
cd demo
test.py
cd ..
cd ..

goto END 

setup.py build --comiler=mingw32 sdist --formats=gztar,zip bdist_wininst


REM Py2exe is not quite automatic. Doubtless this is version and computer specific...
cd py2exe
copy OpenGL.__init__.py c:\python24\lib\site-packages\OpenGl\__init__.py
copy OpenGL.Tk.__init__.py c:\python24\lib\site-packages\OpenGl\Tk\__init__.py
del /q dist build
py2exesetup.py py2exe
copy c:\python24\lib\site-packages\OpenGl\Tk\win32-tk8.4\* dist\
cd dist
del /q _na* _ns* ssl.pyd


:END

