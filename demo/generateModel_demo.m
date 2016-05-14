
%% Create a VariableList object contains the information about the variables
% More detail information can be found in the generateVariableList tutorial
varName = {'pressure','temperature','rate'};
H = repmat([-1,1],3,1);
varVal = zeros(3,1);
varList = generateVar(varName,H,varVal);

%% Create a quadratic model
% To create a B2BDC.B2Bmodels.QModel object by input the coefficient matrix
% and vairable information

% Create a random coefficient matrix
nVar = size(H,1);
coefMatrix = rand(nVar+1);

% Make the coefficient matrix symmetric
coefMatrix = coefMatrix+coefMatrix';

% Create the QModel object
Qmodel = generateModel(coefMatrix,varList);

%% Create a rational quadratic model (1)
% To create a B2BDC.B2Bmodels.RQModel object by input the coefficient matrix
% of both numerator and denominator, vairable information and condition
% number K for the denominator

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

%% Create a rational quadratic model (2)
% To create a B2BDC.B2Bmodels.RQModel object by input the coefficient matrix
% of both numerator and denominator and variable information. The
% positivity condition of the denominator is then checked.

% Create a random coefficient matrix for nominator
Ncoef = rand(nVar+1);

% Make the coefficient matrix symmetric
Ncoef = Ncoef+Ncoef';

% Create a coefficient matrix for denominator that satisfies the positivity
% condition
D1 = zeros(nVar+1);
D1(1,1) = 1;

% Create a coefficient matrix for denominator that does not satisfies the positivity
% condition
D2 = zeros(nVar+1);
D2(1,1) = -1;

% Create the RQModel object with D1
RQmodel = generateModel(Ncoef,D1,varList);

% Create the RQModel object with D2, the following code will fail due to
% the non-positivity of the denominator D2
RQmodel = generateModel(Ncoef,D2,varList);

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

%% Attribution
% UC Berkeley, PSAAP team.  Fall 2015