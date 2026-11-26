function projectDir = setup_project()
% Add all project MATLAB folders and return the repository root.
% From the repository root: addpath('matlab'); setup_project;
projectDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(projectDir,'matlab')));
end
