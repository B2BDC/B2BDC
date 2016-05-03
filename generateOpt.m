function  y = generateOpt(varargin)
% Returns a B2BDC.Option object. The input can be variant number
%  of 'Property_Name' and 'Property_Value' pairs. The valid property
% names and property values are:
% ---------------------------------------------------------------------
   % ConsistencyMeasure:
   %   'relative' -  Calculate consistency measure in the way that
   %                  maximize              y
   %                 subject to  (1-y)*(LB(i)-OB(i)) < M(i)-OB(i)
   %                              M(i)-OB(i) < (1-y)*(UB(i)-OB(i))
   %                 OB(i) is the observation value of ith experiment
   %                 LB(i) is the lower bound of ith experiment
   %                 UB(i) is the upper bound of ith experiment
   %                 M(i) is the surrogate model of ith experiment
   %   'absolute' -  Calculate consistency measure in the way that 
   %                  maximize          y
   %                 subject to   LB(i) + y < M(i)
   %                              M(i) < UB(i) - y
   %   'DClab'    -  Calculate consistency measure from DClab package
% ---------------------------------------------------------------------
   % ExtraLinFraction:
   %   1 ~ 100  -  Include the indicated number of fraction of extra
   %               linear pairs in outer bound calculation
   %      -1    -  Include extra linear pairs with influence factor greater
   %               than 5% of the most influential pair in outer bound
   %               calculation
   % ---------------------------------------------------------------------
   % TolConsis: The tolerance used in bi-sectional algorithm in consistency
   %            mearure outer bound calculation
   % ---------------------------------------------------------------------
   % Display:
   %   true  - Display the progress information in the optimization
   %           procedure
   %   false - Display no progress information in the optimization
   %           procedure
   % ---------------------------------------------------------------------
   % AddFitError:
   %   true  - The calculation in consistency measure and model prediction
   %           will add maximum absolute fitting error of each surrogate
   %           model to the corresponding QOI bound.
   %   false - The calculation in consistency measure and model prediction
   %           will not add maximum absolute fitting error of each surrogate
   %           model to the corresponding QOI bound.
   % ---------------------------------------------------------------------
   % SelfInconsisFlow:
   %    true - When doing the selfInconsisAnalysis, it will show up the
   %           sensitivity plot for every self-inconsistent dataset unit
   %           before showing the final stats
   %   false - Show final stats directly
   % ---------------------------------------------------------------------
   % SOSrelaxOrder:
   %     N   - A positive integer that defines the order of SOS relaxations
   %           in SOS optimization process for polynomial models.
   % ---------------------------------------------------------------------
   % MaxPWsubdom: 
   %     N   - A positive integer that defines the maximum subdomains a
   %           piecewise model can have
   % ---------------------------------------------------------------------
   % PWTol:
   %    tol  - A numerical factor specifies the stop criteria for domain
   %           dividing.


% Created: July 15, 2015    Wenyu Li

if nargin > 0
   y = B2BDC.Option(varargin);
else
   y = B2BDC.Option();
end