function flag = quadratictest(dsUnits)
%   FLAG = QUADRATICTEST(DSUNITS) returns a logical FLAG true if all
%   DSUNITS have surrogate models which are quadratic.
%
%   ********  INTERNAL FUNCTION OF B2BDC ********

% Created June 30, 2015     Wenyu Li

flag = true;
for i = 1:length(dsUnits)
   testModel = dsUnits(i).SurrogateModel;
   if ~isa(testModel,'B2BDC.B2Bmodels.QModel')
      flag = false;
      break;
   end
end
end