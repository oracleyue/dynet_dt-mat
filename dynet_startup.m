% Main Path:
libRootPath = [fileparts(mfilename('fullpath')) '/'];
addpath(libRootPath);


% Additional Paths:

% - for auxiliary scripts
addpath([libRootPath 'auxiliary/']);

% - for simulators
addpath([libRootPath 'simulators/']);
addpath([libRootPath 'simulators/members/'])

% - for the 3rd party solvers
addpath([libRootPath 'solvers/']);
