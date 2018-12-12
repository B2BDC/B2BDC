function y = eval_with_discrepancy(obj,X,DSIdx)
% Y=EVAL_WITH_DISCREPANCY(OBJ,X,DSIDX) evaluates the QOIs, specified by
% DSIdx, or all if DSIdx is not given, with model parameter and discrepancy
% correction coefficients given in X.

%  Created: Sep 12, 2018     Wenyu Li

[idall,Qall] = obj.getQ_RQ_expansion;
if nargin > 2
   if iscell(DSIdx)
      allDS = {obj.DatasetUnits.Values.Name}';
      [~,~,idQOI] = intersect(DSIdx,allDS,'stable');
   else
      [~,~,idQOI] = intersect(DSIdx,1:obj.Length,'stable');
   end
   if length(idQOI) ~= length(DSIdx)
      error('Wrong specified QOI index');
   end
else
   idQOI = 1:obj.Length;
end
nSample = size(X,1);
nQOI = length(idQOI);
y = zeros(nSample,nQOI);
for i = 1:nQOI
   Q0 = Qall{idQOI(i)};
   id0 = idall{idQOI(i)};
   xx = B2BDC.Fitting.expandBasis(X(:,id0));
   vec = B2BDC.Fitting.coef2vec(Q0);
   y(:,i) = xx*vec;
end