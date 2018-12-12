function [QOIrange, QOISensitivity, xOpt, aOpt] = predictQOI_minB(obj, QOIobj, B2Bopt,alpha, QOIcorrection)
% Returns the inner bounds, the feasible points of maximum and minimum 
% values of the QOI model subject to the constraints that both dataset units and
% variables are within their respective bounds and the model discrepancy is
% minimized as |\delta / model|. The sensitivity of the
% outer bounds with respect to experimental bounds and variable bounds is
% also returned.
% The input arguments are:
%    QOIobj - A B2BDC.B2Bmodel.Model object that defines the QOI
%    B2Bopt - B2BDC.Option object (optional)
% The output are:
%    QOIrange - A structure defines the range of QOI subject to the dataset
%    Xopt - A nVar-by-2 double matrix

% Created: Oct 8, 2018   Wenyu Li

if nargin < 3
   B2Bopt = generateOpt('Prediction','inner');
elseif ~isa(B2Bopt,'B2BDC.Option.Option')
   error('Wrong option object')
end
if nargin < 4
   alpha = 0;
end
if nargin > 4 && ~isempty(QOIcorrection)
   q.Model = QOIobj;
   q.Correction = QOIcorrection;
else
   q.Model = QOIobj;
   q.Correction.GroupIndex = 0;
   q.Correction.Value = [];
end
name1 = QOIobj.VarNames;
name2 = obj.VarNames;
[~,id] = intersect(name1,name2);
if length(id) < length(name1)
   rflag = true;
else
   rflag = false;
end
if ~obj.isConsistent(B2Bopt)
   error('The dataset is inconsistent')
else
   tmpCM = obj.ConsistencyMeasure;
   tmpS = obj.ConsistencySensitivity;
   tmpFlag = obj.FeasibleFlag;
   tmpx0 = obj.FeasiblePoint;
end
if rflag
   tVar = QOIobj.Variables;
   ntVar = tVar.Length;
   tCoef = zeros(ntVar+1);
   tCoef(1,1) = 1;
   tM = B2BDC.B2Bmodels.QModel(tCoef,tVar);
   tUnit = generateDSunit('redundant',tM,[0,2]);
   obj.addDSunit(tUnit);
   if ~obj.isConsistent(B2Bopt)
      error('The dataset is inconsistent')
   end
end
polyflag = polytest(obj);
nnflag = networktest(obj);
if polyflag
   if B2Bopt.Display
      disp('=======================================================');
      disp('Calculating QOI range...');
      disp('=======================================================');
   end
   [QOIrange, QOISensitivity, xOpt] = obj.preQOIpoly(QOIobj,B2Bopt,rflag);
   if B2Bopt.Display
      disp('The calculation is done')
      disp(['Minimum value of QOI is within: [' num2str(QOIrange.min(1)) ' ' num2str(QOIrange.min(2)) ']'])
      disp(['Maximum value of QOI is within: [' num2str(QOIrange.max(1)) ' ' num2str(QOIrange.max(2)) ']'])
   end
   return
elseif nnflag
   if B2Bopt.Display
      disp('=======================================================');
      disp('Calculating QOI range...');
      disp('=======================================================');
   end
   [QOIrange, QOISensitivity, xOpt] = obj.preQOInn(QOIobj,B2Bopt,rflag);
   if B2Bopt.Display
      disp('The calculation is done')
      disp(['Minimum value of QOI is within: [' num2str(QOIrange.min(1)) ' ' num2str(QOIrange.min(2)) ']'])
      disp(['Maximum value of QOI is within: [' num2str(QOIrange.max(1)) ' ' num2str(QOIrange.max(2)) ']'])
   end
   return
else
   [minin, maxin, s, xOpt, aOpt] = obj.preQOIfmincon_minB(q,B2Bopt.Display,rflag,B2Bopt,alpha);
   QOIrange.min = minin;
   QOIrange.max = maxin;
%    [minin, maxin, s, xOpt, abE] = obj.preQOIopti(q,B2Bopt.Display,rflag,B2Bopt);
%    toc;
   QOISensitivity.Inner = s;
end

obj.ConsistencyMeasure = tmpCM;
obj.ConsistencySensitivity = tmpS;
obj.FeasibleFlag = tmpFlag;
obj.FeasiblePoint = tmpx0;

end
