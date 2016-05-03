function [varList1, varList2] = subDivide(obj,varIdx,factor)
% Returns two B2BDC.B2Bvariables.VariableList objects. These
% two VariableList objects differ by a single variable range. The
% variable is defined by the varIdx input, it can be either the variable
% name or the variable index in the VariableList. The factor is a real
% value between 0 and 1 and defines the split point xSplit as 
% xSplit = x_LowerBound + factor * (x_UpperBound - x_LowerBound). So the
% variable range is [x_LowerBound, x_Split] in varList1 and [x_Split, x_UpperBound]
% in varList2.
% The syntax of the function is as following:
% [varList1, varList2] = subDivide(pre-varList, varIndex, split_factor)

%  Created: Nov 6, 2015     Wenyu Li

nVar = obj.Length;
if ischar(varIdx)
   [~,id,~] = intersect({obj.Values.Name},varIdx, 'stable');
   if isempty(id)
      error('Input variable is not found')
   end
   if length(id) > 1
      error('Can only delete one variable at one time')
   end
elseif isscalar(varIdx)
   if varIdx > obj.Length
      error('Input variable is not found')
   end
   id = varIdx;
else
   error('Wrong input index type')
end
varList1 = obj;
varList2 = obj;
if factor <= 0 || factor >= 1
   error('Division factor should between 0 and 1')
end
LB = obj.Values(id).LowerBound;
UB = obj.Values(id).UpperBound;
xSplit = LB + factor*(UB-LB);
varList1.Container{id}.UpperBound = xSplit;
varList2.Container{id}.LowerBound = xSplit;