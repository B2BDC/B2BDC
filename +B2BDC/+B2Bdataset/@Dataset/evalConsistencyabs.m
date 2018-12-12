function evalConsistencyabs(obj,b2bopt)
% Calculate the inner and outer bound for the consistency measure (in
% absolute value) of a B2BDC.B2Bdataset.Dataset object. More criteria about the
% optimization process can be manipulated with b2bopt input.

% Created: June 17, 2015      Wenyu Li
% Modified: June 28, 2015      Wenyu Li    (sensitivity added)
% Modified: July 13, 2015     Wenyu Li   (Gradient provided)

frac = b2bopt.ExtraLinFraction;
disflag = b2bopt.Display;
% tic;
[yin,s,xopt,abE,flag] = absCMfmincon(obj,disflag,b2bopt);
% [yin,s,xopt,abE,flag] = absCMopti(obj,disflag,b2bopt);
% toc;
obj.ConsistencySensitivity.Inner = s;

if yin >= 0
   xf = xopt;
   if b2bopt.AddFitError
      obj.FeasibleFlag = true;
   end
%    obj.FeasiblePoint = xf;
elseif disflag
   savexf = input('Do you want to save the optimization point? (y/n) \n','s');
   if strcmpi(savexf,'y')
      save('xopt','xopt');
   end
end

if disflag
   disp('=======================================================');
   disp('Calculating outer bound...');
   disp('=======================================================');
end
if flag
%    [yout,sensitivity] = obj.sedumiconsisquadabs(b2bopt,abE);
   warning('off','all');
   [yout,sensitivity] = obj.cvxconsisquadabs(b2bopt,abE);
   warning('on','all');
   obj.ConsistencyMeasure = [yin yout];
   obj.ConsistencySensitivity.Outer = sensitivity;
else
   %       [yout_result(1),sensitivity{1}] = obj.sedumiconsisabs(yin,b2bopt,abE);
   warning('off','all');
   [yout,sensitivity{1}] = obj.cvxconsisabs(yin,frac,abE);
   warning('on','all');
   obj.ConsistencyMeasure = [yin yout];
   obj.ConsistencySensitivity.Outer = sensitivity;
end
if disflag
   disp(' ')
   disp('The calculation is done')
   disp(['Consistency LB: ' num2str(yin)])
   disp(['Consistency UB: ' num2str(yout)])
end

end