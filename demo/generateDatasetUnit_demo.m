% This is a demo showing how to use generateDSunit function to create
% B2BDC.B2Bdataset.DatasetUnit object.

%  Created: Nov 2, 2015     Wenyu Li

%% Create a string defines the name of the dataset
dsUnitName = 'Example Dataset Unit';

%% Create a model for the dataset
% More detail information can be found in generateModelbyFit demo

% create a VariableList
varName = {'pressure','temperature','rate'};
H = repmat([-1,1],3,1);
varVal = zeros(3,1);
varList = generateVar(varName,H,varVal);

% create a coefficient matrix
nVar = size(H,1);
coefMatrix = rand(nVar+1);
coefMatrix = coefMatrix+coefMatrix';
Qmodel = generateModel(coefMatrix,varList);

%% Create the QoI bounds
QoIBD = [-0.2, 0.45];

%% Create the QoI observed value
QoIval = 0.1;

%% Create a dataset unit
% User can provide a observed value for the QOI or the observed value will
% be set to be the middle point of the bound

% with user input QoI observed value
dsUnit_1 = generateDSunit(dsUnitName,Qmodel,QoIBD,QoIval)

% without user input QoI observed value
dsUnit_2 = generateDSunit(dsUnitName,Qmodel,QoIBD)

%% Attribution
% UC Berkeley, PSAAP team.  Fall 2015