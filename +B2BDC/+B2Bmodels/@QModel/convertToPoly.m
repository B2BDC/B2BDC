function pModel = convertToPoly(obj)
% PMODEL = CONVERTTOPOLY(OBJ) converts a QModel object into a PolyModel
% object.

%  Created: Dec 16, 2016     Wenyu Li

m1 = obj.CoefMatrix;
c1 = B2BDC.Fitting.coef2vec(m1);
vars = obj.Variables;
s1 = generateSupport2(length(c1),vars.Length);
er = obj.ErrorStats;
pModel = B2BDC.B2Bmodels.PolyModel(s1,c1,vars,er);



   function y = generateSupport2(nS,nV)
      y = zeros(nS,nV);
      y(2:nV+1,:) = eye(nV);
      tM = eye(nV);
      tM(:,1) = tM(:,1)+1;
      count = nV+2;
      for i = 1:nV
         y(count:count+nV-i,i:end) = tM;
         tM(:,end) = [];
         tM(end,:) = [];
         count = count + nV + 1 - i;
      end
   end


end