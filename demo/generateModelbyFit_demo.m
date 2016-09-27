% This is a demo showing how to use generateModelbyFit function to create
% B2BDC.B2Bmodels.Model object with different input formats.

%  Created: Nov 1, 2015     Wenyu Li

%% Create a VariableList object contains the information about the variables
% More detail information can be found in the generateVariable demo

varName = {'pressure','temperature','rate'};
H = repmat([-1,1],3,1);
varList = generateVar(varName,H);

%% Create X data by sampling in the variable domain
% The X data should be a nSample-by-nVar matrix

% number of variables
nVar = size(H,1);

% number of sample points
nSample = 20*nVar;

% generate LHS sample points from the VariableList object
xData = varList.makeLHSsample(nSample);

%% Create Y data from the function handle
% Create a function handle
f = @(x)exp(x(:,1).*x(:,2))+x(:,1).^2 - x(:,1).*x(:,2)+3*x(:,2);

% generate nSample-by-1 column vector Y data
yData = f(xData);

%% Create a B2BDC.B2Bmodels.QModel by fitting the data
% Two different criteria to minimize different norm of error terms can be used
% by a string input shown in the following example:

% generate a linear model by minimize 2-norm error
Lmodel = generateModelbyFit(xData,yData,varList,'lin')

% generate a quadratic model by minimize infinity-norm error
Qmodel_1 = generateModelbyFit(xData,yData,varList,'qinf')

% generate a quadratic model by minimize 2-norm error
Qmodel_2 = generateModelbyFit(xData,yData,varList,'q2norm')

%% Create a B2BDC.B2Bmodels.RQModel by fitting the data
% To generate a RQmodel by fitting the data, user can provide a positive
% condition number K, or the algorithm will find a K between 1-20 by cross
% validation

% generate the option object
opt = generateOpt('Display',false);

% The condition number is calculated by cross-validation
RQmodel_1 = generateModelbyFit(xData,yData,varList,'rq',opt)

% The condition number is given by user
RQmodel_2 = generateModelbyFit(xData,yData,varList,'rq',opt,10)

%% Attribution
% UC Berkeley, PSAAP team.  Fall 2015