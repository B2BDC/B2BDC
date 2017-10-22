% This is a demo for how to create a B2BDC.B2Bmodels.Model object by
% direct input and some function examples of the object.More detail 
% information about inputs of the function 'generateQModel' and 
% 'generateRQModel' can be found by typing help generateQModel and 
% generateRQMode in command window

%  Created: Nov 2, 2015     Wenyu Li

%% Create a VariableList object contains the information about the variables
% More detail information can be found in the VariableList_demo
varName = {'pressure','temperature','rate'};
H = repmat([-1,1],3,1);
varVal = zeros(3,1);
varList = generateVar(varName,H,varVal);

%% Create a quadratic model
% Create a random coefficient matrix
nVar = size(H,1);
coefMatrix = rand(nVar+1);

% Make the coefficient matrix symmetric
coefMatrix = coefMatrix+coefMatrix';

% Create the QModel object
Qmodel = generateModel(coefMatrix,varList);

%% Create a rational quadratic model
% Create a random coefficient matrix for nominator
Ncoef = rand(nVar+1);

% Make the coefficient matrix symmetric
Ncoef = Ncoef+Ncoef';

% Create a random coefficient matrix for denominator
Dcoef = rand(nVar+1);

% Make the coefficient matrix symmetric
Dcoef = Dcoef+Dcoef';

% Set the denominator outer bound K value
k = 10;

% Create the RQModel object
RQmodel = generateModel(Ncoef,Dcoef,varList,k);

%% Evaluate a function value at given x values
% The input should be a nSample-by-nVar matrix that defines the values of x
% you want to calculate the function value. The output is a nSample-by-1
% column vector of the function values at those points.

% Make sample points, more information can be found in VariableList_demo
nSample = 20;
xSample = varList.makeLHSsample(nSample);

% Calculate the function values from the QModel object
y_q = Qmodel.eval(xSample);

% Calculate the function values from the RQModel object
y_rq = RQmodel.eval(xSample);