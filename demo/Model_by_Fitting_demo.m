% This is a demo for how to create a B2BDC.B2Bmodels.Model object by
% fitting data. More detail information about inputs of the function
% 'generateModelbyFit' can be found by typing help generateModelbyFit in
% command window

%  Created: Nov 1, 2015     Wenyu Li

%% Create a VariableList object contains the information about the variables
% More detail information can be found in the VariableList_demo
varName = {'pressure','temperature','rate'};
H = repmat([-1,1],3,1);
varVal = zeros(3,1);
varList = generateVar(varName,H,varVal);

%% Create a random X data
nVar = size(H,1);
nSample = 20*nVar;
xData = varList.makeLHSsample(nSample);

%% Create a random Y data
yData = rand(nSample,1);

%% Create a B2BDC.B2Bmodels.QModel by fitting the data
% generate an option object
b2bopt = generateOpt('Display',false);

% generate a quadratic model by minimize infinity-norm error
Qmodel_1 = generateModelbyFit(xData,yData,varList,'qinf',b2bopt);

% generate a quadratic model by minimize 2-norm error
Qmodel_2 = generateModelbyFit(xData,yData,varList,'q2norm',b2bopt);

%% Create a B2BDC.B2Bmodels.QModel by fitting the data
% generate a rational quadratic model by minimize infinity-norm error

% The upper bound of denominator is calculated by cross-validation
RQmodel_1 = generateModelbyFit(xData,yData,varList,'rq',b2bopt);

% The upper bound of denominator is given by user
RQmodel_2 = generateModelbyFit(xData,yData,varList,'rq',b2bopt,10);