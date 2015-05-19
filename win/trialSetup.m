%Add Paths

%Get absolute path to the folder containing this file
basePath = [fileparts(mfilename('fullpath')) filesep];


%GAMPMATLAB paths
addpath([basePath '../../main']) %main GAMPMATLAB code
addpath([basePath './comparison/l1magic/Optimization']);
addpath([basePath './comparison/ApprocMsgPassing']);
addpath([basePath './comparison/sparseLab']);

%Handle random seed
defaultStream = RandStream.getGlobalStream;
if 1
    savedState = defaultStream.State;
    save random_state.mat savedState;
else
    load random_state.mat %#ok<UNRCH>
end
defaultStream.State = savedState;



