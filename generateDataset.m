function y = generateDataset(dsName)
% Returns a B2BDC.B2Bdataset.Dataset object with no dataset units.
% The input argument is:    
%
%   dsName - A string defines the name of the dataset

%  Created: Nov 1, 2015    Wenyu Li

if nargin > 1
   error('Wrong number of inputs');
end
if nargin == 1 && ischar(dsName)
   y = B2BDC.B2Bdataset.Dataset(dsName);
else
   y = B2BDC.B2Bdataset.Dataset('Dataset 1');
end