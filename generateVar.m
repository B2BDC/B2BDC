function y = generateVar (varName, H, val)
%  Returns a VariableList object from array of Names and Bounds
%
% The input arguments are:
% 
%  varName - cell array of variable names
%     H    - A n_variable-by-2 matrix defines the domain of the variable
%    val   - Optional, A n_variable-by-1 column vector specifies the nominal 
%            value of the variables. If not input by user, will calculate
%            as the mean of H for each variable
%
% The input formats are:
%    V = generateVar(NAMES,BOUNDS) creates a VariableList object from an
%    N-by-1 cell array NAMES which defines the name of each variable, and 
%    N-by-2 numeric array BOUNDS which defines the interval for each variable.  
%    The nominal value of each parameter is defined as the midpoint of the
%    interval specified in BOUNDS
%
%    V = generateVar(NAMES,BOUNDS,NOMINAL) sets the nominal value for each
%    parameter as the corresponding element in the N-by-1 array NOMINAL.

% Created: July 13, 2015     Wenyu Li
%  Modified: Jan 13, 2016     Wenyu Li (comments edited)

n_variable = length(varName);
if n_variable ~= size(H,1)
   error('generateVar:InputAlignment',...
       'The input variable name and variable bounds have different dimension.')
end
y = B2BDC.B2Bvariables.VariableList();
for i = 1:n_variable
    variableName = varName{i};
    LB = H(i,1);
    UB = H(i,2);
    if UB <= LB
       error('generateVar:IllegalBounds',...
           ['Invalid lower/upper bounds for ' variableName]);
    end
    if nargin > 2
        nominalVal = val(i);
    else
        nominalVal = 0.5*(LB+UB);
    end
    modelVar = B2BDC.B2Bvariables.ModelVariable(variableName,LB,UB,nominalVal);
    y = y.add(modelVar);
end
