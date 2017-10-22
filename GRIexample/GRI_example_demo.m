% B2BDC example: GRI-Mech 3.0

% Created: June 26, 2015    myf, Wenyu Li
% Modified: December 8, 2015 Team

%% Input: experimental QOI data (user-supplied)
% The experimental QOI data is a user-supplied Excel file of the form, exemplified below, 
%
% <<..\Images\ExcelSnapshot.png>>
% 
% where the columns, in order, are
%
% * A: |name| - unique name identifying experimental QOI (e.g., 'bch2o.t1')       
% * B: |LB|  - LOWER bound on observed value of QOI                          
% * C: |UB|  - UPPER bound on observed value of QOI                          
% * D: |value| - observed value of QOI (optional)                            
%
[~,~,experimentData] = xlsread('GRIMech_expdata.xls');
% creating Matlab cell array (numQOI-by-numCol) from the Excel file.
% 
%% Input: 
%create input data (cell array) of models "reproducing" QOIs

% modelData = { name M vars } is a numberOfObservations-by-3 cell array
%    key - unique key identifying experimental QoI (same as in the experimentalData)
%    M   - surrogate model specified by coefficients or by data to fit (see below)
%    vars - model variables, a numberOfVariables-by-4 cell array as { varName LB UB x0 }
%         varname - variable key (e.g., 'k1' or 'A(O+H2=H+OH)')
%         LB  - LOWER bound of the variable
%         UB  - UPPER bound of the variable
%         x0  - nominal value of the variable (typically x0 = (LB + UB)/2)
% ________________________________________________________________________
%
% (I) M is specified by model coefficients
%
%  (a) Surrogate model is a quadratic model
%      M is a numberOfObservations-by-2 cell array, M = { key coef }
%      key = 'quadratic'
%      coef = [    c0   0.5*c1   0.5*c2  ...
%              0.5*c1       c11  0.5*c12 ...
%                ...                         ]
%           - it is a quadratic form that evaluates as [1;x]' * coef * [1;x]
%             where x is a column-vector of variables
%
%  (b) Surrogate model is a rational-quadratic model
%      M is a numberOfObservations-by-4 cell array, M = { key numeratorCoef denominatorCoef kVal }
%      key is 'rational quadratic'
%      numeratorCoef and denominatorCoef are quadratic forms, similar to coef of (a)
%         such that M evaluates to
%            [1;x]' * numeratorCoef * [1;x] / [1;x]' * denominatorCoef * [1;x]
%      kVal is a scalar greater than 1 that it is an upper bound on the value of the denominator
%            1 <= [1;x]' * denominatorCoef * [1;x] <= kVal
%
% (II) M is specified by data to fit
%      M is a numberOfObservations-by-3 cell array, M = { key, Xdata, Ydata }
%      key is the same as in (I) above
%      Xdata and Ydata are related by Ydata = fQOI(Xdata), where
%         fQOI is a function that reproduces the corresponding QoI
%         Xdata is the input (design) matrix to fQOI; 
%               it is a nSample-by-numberOfVariables double array, whose
%               columns are variables specified by 'vars' in modelData 
%               and rows specify individual runs of fQOI
%         Ydata is a nSample-by-1 double array of computed QoI values

s = load('GRIMech_modeldata.mat');
modelData = s.GRIMech_modeldata;

%% B2B-DC begins: Create the GRI Mech 3.0 dataset
dsGRI = generateDataset('GRI Mech 3.0');

%% Add dataset units into the dataset
dsGRI = addData(dsGRI,experimentData,modelData);
