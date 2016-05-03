function  y = eval(obj,x, varObj)
%   Y = EVAL(OBJ, X) evaluates a surrogate model of a DatasetUnit at X
%   sampled points to produce a column vector Y of the surrogate model's
%   output. X is a matrix of size nSample-by-nVariable, where nSample is
%   the number of samples to be evaluated and nVariable is the number of
%   variables in the dataset unit.
%
%   Y = EVAL(OBJ, X, VAROBJ) also evaluates a surrogate model of a
%   DatasetUnit at X sampled points, where X can be of size greater than
%   nVariable. VAROBJ will specify which columns of X to be used in
%   evaluating the surrogate model to produce a column vector Y.
%
%   ********  INTERNAL FUNCTION OF B2BDC ********

% Created: July 21, 2015    Wenyu Li

unitModel = obj.SurrogateModel;
if nargin > 2
   y = unitModel.eval(x, varObj);
else
   y = unitModel.eval(x);
end
