@echo on
@REM mkdir build
cd build
cmake --build . --config Debug
cd .. 

cd dll
xcopy "*" "..\build\Debug\" /Y
cd ..