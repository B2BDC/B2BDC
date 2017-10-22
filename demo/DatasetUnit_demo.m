% This is a demo for how to create a B2BDC.B2Bdataset.DatasetUnit
% object. More detail information can be found by typing help
% generateDSunit in command window

%  Created: Nov 2, 2015     Wenyu Li

%% Create a string defines the name of the dataset
dsUnitName = 'Example Dataset Unit';

%% Create a model for the dataset
% More detail information can be found in Model_by_Input_demo

% Create a VariableList
varName = {'pressure','temperature','rate'};
H = repmat([-1,1],3,1);
varVal = zeros(3,1);
varList = generateVar(varName,H,varVal);

% Create a coefficient matrix
nVar = size(H,1);
coefMatrix = rand(nVar+1);
coefMatrix = coefMatrix+coefMatrix';
Qmodel = generateModel(coefMatrix,varList);

%% Create the QoI bounds
QoIBD = [-0.2, 0.45];

%% Create the QoI observed value
QoIval = 0.1;

%% Create a dataset unit
% With user input QoI observed value
dsUnit_1 = generateDSunit(dsUnitName,Qmodel,QoIBD,QoIval);

% Without user input QoI observed value
dsUnit_2 = generateDSunit(dsUnitName,Qmodel,QoIBD);