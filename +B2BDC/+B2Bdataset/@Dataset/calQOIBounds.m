function qoiPos = calQOIBounds(obj, index, opt)
%   QOIPOS = CALQOIBOUNDS(OBJ) returns a structure containing 
%   the calculated bounds of all DatasetUnits inside the dataset OBJ,
%   where nDSUnit is the number of DatasetUnits inside the dataset.
%
%   QOIPOS = CALQOIBOUNDS(OBJ, INDEX) returns a structure containing 
%   the bounds of DatasetUnits specified by INDEX. INDEX is either a 
%   cell array of DatasetUnit names or a vector of indicies of the n number
%   of DatasetsUnits to be calculated. 
%   
%   QOIPOS = CALQOIBOUNDS(OBJ, INDEX, OPT) the OPT input can be used to 
%   modify the printed output or optimization criteria. OPT is a 
%   B2BDC.Option object.

%  Created: Nov 30, 2015     Wenyu Li
%  Modified: Nov 16, 2016     Wenyu Li

if nargin < 2
   index = [];
end
if nargin < 3
   opt = generateOpt;
   opt.Display = false;
   opt.ExtraLinFraction = -1;
end
pflag = opt.Prediction;
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
selectName = allName(id-1);
qbd1 = zeros(length(id),2);
qbd2 = zeros(length(id),2);
% if opt.Display
   h = waitbar(0,'Calculating the Posterior QOI Bounds');
% end
for i = 1:length(id)
   tep_QOI = obj.DatasetUnits.Values(id(i)).SurrogateModel;
   QOIrange = obj.predictQOI(tep_QOI, opt);
   if strcmp(pflag,'both')
      qbd1(i,1) = QOIrange.min(1);
      qbd1(i,2) = QOIrange.max(2);
      qbd2(i,1) = QOIrange.min(2);
      qbd2(i,2) = QOIrange.max(1);
   else
      qbd1(i,1) = QOIrange.min;
      qbd1(i,2) = QOIrange.max;
   end
%    if opt.Display
      waitbar(i/length(id),h);
%    end
end
% if opt.Display
   close(h);
% end
qoiPos.varName = selectName;
switch pflag
   case 'both'
      qoiPos.OuterBound = qbd1;
      qoiPos.InnerBound = qbd2;
   case 'inner'
      qoiPos.InnerBound = qbd1;
   case 'outer'
      qoiPos.OuterBound = qbd1;
end