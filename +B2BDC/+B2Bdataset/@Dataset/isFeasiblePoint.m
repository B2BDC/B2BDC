function isFeas = isFeasiblePoint(ds, point)
%   ISFEAS = ISFEASIBLEPOINT(DS, POINT, OPT) returns a logical true if POINT is 
%   feasible for the dataset DS. POINT is a nPoint-by-nVariables vector of 
%   the points which you want to examine feasibility. If any points are
%   infeasible, a logical false is returned.

nVar = ds.Variables.Length;
if size(point,2) ~= nVar
   point = point';
end
if size(point,2) ~= nVar
   error('Invalid input sample point size')
end
nPoint = size(point,1);
isFeas = true(nPoint,1);

dsBounds = ds.calBound;
if ds.FeasibleFlag
   for i = 1:ds.Length
      abE = ds.DatasetUnits.Values(i).SurrogateModel.ErrorStats.absMax;
      dsBounds(i,1) = dsBounds(i,1) - abE;
      dsBounds(i,2) = dsBounds(i,2) + abE;
   end
end
% ic = 0;
nmax = 10^5;
n1 = floor(nPoint/nmax);
n2 = mod(nPoint,nmax);

for i1 = 1:n1
   dsEval = ds.eval(point(1:nmax,:));
%    for i = 1:nmax
%       ic = ic+1;
%       if any(dsEval(i,:) < dsBounds(:,1)')
%          isFeas(ic) = false;
%          continue
%       elseif any(dsEval(i,:) > dsBounds(:,2)')
%          isFeas(ic) = false;
%          continue
%       elseif ~ds.Variables.isFeasiblePoint(point(i,:))
%          isFeas(ic) = false;
%          continue
%       end
%    end
   a1 = dsEval > repmat(dsBounds(:,1)',nmax,1);
   a2 = dsEval < repmat(dsBounds(:,2)',nmax,1);
   b1 = all(a1');
   b2 = all(a2');
   b3 = ds.Variables.isFeasiblePoint(point(1:nmax,:));
   id = b1 & b2 & b3';
   id = [true(1,(i1-1)*nmax), id];
   isFeas(~id) = false;
   point(1:nmax,:) = [];  
end

if n2 > 0
dsEval = ds.eval(point);

% for i = 1:n2
%    ic = ic+1;
%    if any(dsEval(i,:) < dsBounds(:,1)')
%       isFeas(ic) = false;
%       continue
%    elseif any(dsEval(i,:) > dsBounds(:,2)')
%       isFeas(ic) = false;
%       continue
%    elseif ~ds.Variables.isFeasiblePoint(point(i,:))
%       isFeas(ic) = false;
%       continue
%    end
a1 = dsEval >= repmat(dsBounds(:,1)',n2,1);
a2 = dsEval <= repmat(dsBounds(:,2)',n2,1);
b1 = all(a1');
b2 = all(a2');
b3 = ds.Variables.isFeasiblePoint(point(1:n2,:));
id = b1 & b2 & b3';
id = [true(1,n1*nmax), id];
isFeas(~id) = false;
% end
end
