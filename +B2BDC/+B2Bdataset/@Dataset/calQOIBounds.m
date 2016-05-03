function qoiBound = calQOIBounds(obj, index, opt)
%   QOIBOUND = CALQOIBOUNDS(OBJ) returns a nDSUnit-by-2 matrix containing 
%   the calculated outer bounds of all DatasetUnits inside the dataset OBJ,
%   where nDSUnit is the number of DatasetUnits inside the dataset. The 
%   outer bounds are constrained such that the dataset units and 
%   variables are within their corresponding bounds. 
%
%   QOIBOUND = CALQOIBOUNDS(OBJ, INDEX) returns a n-by-2 matrix containing 
%   the outer bounds of DatasetUnits specified by INDEX. INDEX is either a 
%   cell array of DatasetUnit names or a vector of indicies of the n number
%   of DatasetsUnits to be calculated. 
%   
%   QOIBOUND = CALQOIBOUNDS(OBJ, INDEX, OPT) the OPT input can be used to 
%   modify the printed output or optimization criteria. OPT is a 
%   B2BDC.Option object.

%  Created: Nov 30, 2015     Wenyu Li

if nargin < 2
   index = [];
end
if nargin < 3
   opt = generateOpt;
   opt.Display = false;
end
nUnit = obj.Length;
if index > nUnit
    error('Index cannot exceed the number of QOIs in the dataset.');
end
allName = {obj.DatasetUnits.Values.Name};
if isempty(index)
   id = 1:nUnit;
elseif iscell(index)
   [~,~,id] = intersect(index, allName, 'stable');
elseif isnumeric(index)
   id = index;
else
   error('Wrong input index format')
end
qoiBound = zeros(length(id),2);
% if opt.Display
   h = waitbar(0,'Calculating the Posterior QOI Bounds');
% end
for i = 1:length(id)
   tep_QOI = obj.DatasetUnits.Values(id(i)).SurrogateModel;
   QOIrange = obj.predictQOI(tep_QOI, opt);
   qoiBound(i,1) = QOIrange.min(1);
   qoiBound(i,2) = QOIrange.max(2);
%    if opt.Display
      waitbar(i/length(id),h);
%    end
end
% if opt.Display
   close(h);
% end