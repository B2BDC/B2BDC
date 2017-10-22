% This is a demo showing how to create a B2BDC.B2Bvariables.VariableList
% object and its commonly used methods.

%  Created: Nov 1, 2015     Wenyu Li

%% Create a cell array containing variable names
varNames = {'pressure','temperature','rate'};
%% Create a numberOfVaribles-by-2 matrix containing the lower and upper 
% bounds of the variables
varBounds = repmat([-1,1],[3,1]);
%% Create a numberOfVariables-by-1 column vector containing the nominal 
% values of the variables
varNominalVals = zeros(3,1);
%% Create a |VariableList| object
varList1 = generateVar(varNames,varBounds,varNominalVals);
%% Verify the class of varList1
class(varList1)
%% Verify varList1 consists of three variables
length(varList1)
%% Append a VariableList to another VariableList
% The |addList| method can append a |VariableList| object to another
% |VariableList|
% First let's create another |VariableList| object
varNames = {'time','length'};
varBounds = repmat([-2,2],[2,1]);

% If no nominal value input is provided in |generateVar|, we will
% automatically calculate the mean of variable lower and upper bounds 
% as the nominal value
varList2 = generateVar(varNames,varBounds);

% Append varList2 to varList1
varList3 = varList1.addList(varList2);

% Verify varList3 consists of five variables
varList3.Length

%% Create a sample by latin-hypercube design in the variable space
% This function creates a nSample-by-nVar matrix that samples the
% parameter space by latin-hypercube-design method
nSample = 20;
varSample = varList3.makeLHSsample(nSample);

% varSample is 20-by-5
size(varSample)

%% Conclusion
% In this demo, we saw how to create a |VariableList| object and how to 
% appened another |VariableList| object together.  We also saw how to 
% sample the domain defined by the |VariableList| object. 
% More information can be found with |help|. 
help generateVar

%% Attribution
% UC Berkeley, PSAAP team.  Fall 2015