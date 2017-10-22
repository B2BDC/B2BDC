function  y = generatePOPOpt(varargin)
% Returns a B2BDC.Option.POPOption object. The input can be variant number
%  of 'Property_Name' and 'Property_Value' pairs. The valid property
% names and property values are:
   % ---------------------------------------------------------------------
   % relaxOrder:
   %   A integer indicates the relaxation order of the method
   % ---------------------------------------------------------------------
   % sparseSW:
   %   0 - Dense SOS relaxation
   %   1 - Sparse SOS relaxation
   % ---------------------------------------------------------------------
   % perturbation:
   %   A small positive number used to overcome multiple optimal situation
   %   of the original POP
   % ---------------------------------------------------------------------
   % SDPsolver:
   %   sedumi
   %   sdpa
   % ---------------------------------------------------------------------
   % printFileName:
   %   0 - no solution information display
   %   1 - display solution information
   %  'filename' - print out solution information
   % ---------------------------------------------------------------------
   % printLevel:
   % ---------------------------------------------------------------------
   % POPsolver:
   %   active-set
   %   trust-region-reflective
   %   interior-point
   % ---------------------------------------------------------------------
   % POPsolver:
   %   0 - not emplementing mex file
   %   1 - emplementing mex file

% Created: Dec 15, 2016    Wenyu Li

if nargin > 0
   y = B2BDC.Option.POPOption(varargin);
else
   y = B2BDC.Option.POPOption();
end