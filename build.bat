@echo off

set CommonCompilerFlags= -MT -nologo -GR- -Gm- -EHa- -Oi -Od -WX -W4 -wd4100 -wd4127 -wd4101 -wd4201 -wd4189 -wd4505 -DDEBUG_MODE=1 -DSSE=0 -FC -Z7 

mkdir "build"
pushd "build"
cl %CommonCompilerFlags% -O2 -DSSE=0 -DSHADING=0 -DGRAVITY=0 ..\win_particle.cpp User32.lib Gdi32.lib
