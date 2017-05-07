% This is a demo to show how to create a B2BDC.B2Bvariables.VariableList
% object and examples of methods this object has.

%  Created: Nov 1, 2015     Wenyu Li

%% Create a B2BDC.B2Bvariables.VariableList object
% More information can be found by typing help generateOpt in the command
% window

% Create a cell array containing the variable names
varName = {'pressure','temperature','rate'};

% Create a nVar-by-2 matrix containing the lower and upper bounds of the
% variables
H = repmat([-1,1],3,1);

% Create a nVar-by-1 column vector containing the nominal values of the
% variables
varVal = zeros(3,1);

% Create a VariableList object
varList1 = generateVar(varName,H,varVal);

%% Add a VariableList to another VariableList
% This function can merge(combine) a VariableList object to another one

% Create another VariableList object
varName = {'time','length'};
H = repmat([-2,2],2,1);
% No nominal value input, it will calculate the mean of variable lower and
% upper bounds as the nominal value

varList2 = generateVar(varName,H);

% Add varList2 to varList1
varList3 = varList1.addList(varList2);

%% Create a sample by latin-hypercube design in the variable space
% This function can create a nSample-by-nVar matrix that samples the
% parameter space by latin-hypercube-design method
n_Sample = 20;
x_Sample = varList3.makeLHSsample(n_Sample);