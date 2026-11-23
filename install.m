function install

cocodir = fileparts(mfilename('fullpath'));
addpath(fullfile(cocodir, 'toolbox'));
addpath(fullfile(cocodir, 'toolbox','misc'));
addpath(fullfile(cocodir, 'toolbox','MTSHT_ODEsolvers'));

GEFinstall = fullfile(cocodir, 'ext', 'GEF_quat','startup.m');
run(GEFinstall);

end