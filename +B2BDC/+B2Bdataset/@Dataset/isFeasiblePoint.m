function isFeas = isFeasiblePoint(ds, point)
%   ISFEAS = ISFEASIBLEPOINT(DS, POINT) returns a logical true if POINT is 
%   feasible for the dataset DS. POINT is a nPoint-by-nVariables vector of 
%   the points which you want to examine feasibility. If any points are
%   infeasible, a logical false is returned.

nPoint = size(point,1);
isFeas = true(nPoint,1);

dsBounds = ds.calBound;
xBound = ds.Variables.calBound;

dsEval = ds.eval(point);

for i = 1:nPoint
   if any(dsEval(i,:) < dsBounds(:,1)')
      isFeas(i) = false;
   elseif any(dsEval(i,:) > dsBounds(:,2)')
      isFeas(i) = false;
   elseif any(point(i,:) < xBound(:,1)')
      isFeas(i) = false;
   elseif any(point(i,:) > xBound(:,2)')
      isFeas(i) = false;
   end
end

