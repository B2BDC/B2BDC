function  y = generateOptimOpt(varargin)
% Returns a B2BDC.Option.OptimOption object. The input can be variant number
%  of 'Property_Name' and 'Property_Value' pairs. The valid property
% names and property values are:
   % ---------------------------------------------------------------------
   % OptimizationMethod:
   %   LSF(default) - Least square optimization in feasible set
   %   LSH - Least square optimization in prior H
   %   1NF - 1 norm optimization in F
   % ---------------------------------------------------------------------
   % PenaltyWeight:
   %   relative(default) - squared relative difference
   %   absolute - squared absolute difference
   %   user-defined - user defined weights
   % ---------------------------------------------------------------------


% Created: Oct 11, 2016    Wenyu Li

if nargin > 0
   y = B2BDC.Option.OptimOption(varargin);
else
   y = B2BDC.Option.OptimOption();
end