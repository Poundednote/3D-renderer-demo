@echo off

set CommonCompilerFlags= -MT -nologo -GR- -Gm- -EHa- -Oi -Od -WX -W4 -wd4100 -wd4101 -wd4201 -wd4189 -wd4505 -DDEBUG_MODE=1 -FC -Z7 

mkdir "w:\C++\PhysicsEngine\build"
pushd "w:\C++\PhysicsEngine\build"
cl %CommonCompilerFlags% ..\win_particle.cpp User32.lib Gdi32.lib
popd
