function evalConsistencyabsPoly(obj,b2bopt)
% Calculate the inner and outer bound for the consistency measure (in
% absolute value) of a B2BDC.B2Bdataset.Dataset object in the situation that
% all surrogate models are polynomial models and at least one of them is not
% quadratic. More criteria about the optimization process can be manipulated 
% with b2bopt input.

%  Created: Dec 22, 2016     Wenyu Li

disflag = b2bopt.Display;
[yin,s,xopt,abE] = absCMfminconPoly(obj,disflag,b2bopt);
obj.ConsistencySensitivity.Inner = s;

if yin >= 0
   xf = xopt;
   if b2bopt.AddFitError
      obj.FeasibleFlag = true;
   end
   obj.FeasiblePoint = xf;
end

if disflag
   disp('=======================================================');
   disp('Calculating outer bound...');
   disp('=======================================================');
end

[yout,sensitivity] = obj.sosconsispolyabs(b2bopt,abE,b2bopt.SOSOption);
obj.ConsistencyMeasure = [yin yout];
obj.ConsistencySensitivity.Outer = sensitivity;