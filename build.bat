@echo off
call E:\VisualStuio\Community\VC\Auxiliary\Build\vcvarsall x64
set CommonCompilerFlags= -MT -nologo -GR- -EHa- -Oi -Od -wd4100 -wd4127 -wd4101 -wd4201 -wd4189 -wd4505 -w -DDEBUG_MODE=1 -DSSE=0 -FC -Z7
mkdir "build"
pushd "build"
cl %CommonCompilerFlags% -O2 -DPHYSICS_TEST=0 -DPOST_FX=1 -DSSE=0 -DSHADING=1 -DGRAVITY=0 ..\win_particle.cpp User32.lib Gdi32.lib
