function networkConsisabs(obj,b2bopt)
% Calculate the inner bound for the consistency measure (in 
% percentage) of a B2BDC.B2Bdataset.Dataset object. More criteria about the
% optimization process can be specified with b2bopt input.

% Created: April 27, 2018      Wenyu Li


disflag = b2bopt.Display;
[yin,s,xopt] = absCMfminconNN(obj,disflag,b2bopt);
obj.ConsistencySensitivity.Inner = s;
obj.ConsistencySensitivity.Outer = [];
if yin >= 0
   xf = xopt;
   if b2bopt.AddFitError
      obj.FeasibleFlag = true;
   end
   obj.FeasiblePoint = xf;
elseif disflag
   savexf = input('Do you want to save the optimization point? (y/n) \n','s');
   if strcmpi(savexf,'y')
      save('xopt','xopt');
   end
end
yout = inf;
obj.ConsistencyMeasure = [yin yout];
if disflag
   disp(' ')
   disp('The calculation is done')
   disp(['Consistency LB: ' num2str(yin)])
   disp(['Consistency UB: ' num2str(yout)])
end