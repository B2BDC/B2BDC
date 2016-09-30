function install_b2bdc()

% Bound-to-Bound Data Collaboration (B2BDC) installation script
%
%

w = what('+B2BDC');
d = fileparts(w.path);
addpath(fullfile(d));
savepath
