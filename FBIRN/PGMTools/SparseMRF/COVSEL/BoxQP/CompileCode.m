% Script to compile mex file using MATLAB interface

% Test for mex extension to determine matlab version and platform
switch mexext
    case {'mexw32'} % Win32 MATLAB after 7.1
        files='BoxQP.c BoxQP_mex.c utils.c';
        libs=[matlabroot,'\extern\lib\win32\microsoft\libmwlapack.lib'];
        switches='-v -O -DWIN32';
        eval(['mex ',switches,' ',files,' ','''',libs,'''']);
    case {'dll'} % Win32 MATLAB before 7.1
        files='BoxQP.c BoxQP_mex.c utils.c';
        libs=[matlabroot,'\extern\lib\win32\microsoft\msvc60\libmwlapack.lib'];
        switches='-v -O -DWIN32';
        eval(['mex ',switches,' ',files,' ','''',libs,'''']);
    case {'mexmac'}% Macintosh using the VecLib framework
        files='BoxQP.c BoxQP_mex.c utils.c';
        switches='-v -O -Dmac';
        eval(['mex ',switches,' ',files,' ']);
    otherwise % All other platforms
        files='BoxQP.c BoxQP_mex.c utils.c';
        switches='-v -O -Dlinuxp';
        eval(['mex ',switches,' ',files,' ']);
end
disp(' ......................... Done Compiling Source .........................')