function  y = generateSampleOpt(varargin)
% Returns a B2BDC.Option.SampleOption object. The input can be variant number
%  of 'Property_Name' and 'Property_Value' pairs. The valid property
% names and property values are:
   % ---------------------------------------------------------------------
   % SampleMethod:
   %   Sobol - Rejection sampling with sobolset on rotated box
   %   HR - Rejection sampling with hit and run in general polytope
   % ---------------------------------------------------------------------
   % SampleApproximation:
   %   none - No apporximation method is used
   %   DR - Dimension reduction approximation
   %   UA - Uncertainty approximation
   %   DR&UA - Dimension reduction and uncertainty approximation
   % ---------------------------------------------------------------------
   % BatchMaxSample:
   %   The maximum samples collected in one batch
   % ---------------------------------------------------------------------
   % StepInterval:
   %   The step interval in hit and run algorithm
   % ---------------------------------------------------------------------
   % PCTruncation:
   %   The threshold variance value for principal direction truncation
   % ---------------------------------------------------------------------
   % ExtraCut:
   %   The number of extra random direction cut used to generate containing
   %   polytope


% Created: Oct 11, 2016    Wenyu Li

if nargin > 0
   y = B2BDC.Option.SampleOption(varargin);
else
   y = B2BDC.Option.SampleOption();
end