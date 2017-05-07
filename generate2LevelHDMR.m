function pModel = generate2LevelHDMR(f, vars, d, x0)
% Returns a B2BDC.Models.PolyModel object following the cutting HDMR criteria. 
% Only 2nd order variable interaction is considered with maximum d degree
% monomials.
%  f - Function handle that takes in only one input design matrix of size
%      nSample-by-nVar. The output is a column vector of length nSample.
%  vars - A VariableList object specifies the variable information for f
%  d - The maximum degree of polynomials used in the fitting
%  x0 - Optional input cutting point, if not specified, the mean is used

%  Created: Jan 6, 2016     Wenyu Li

n1 = 15;
nVar = vars.Length;
np = 1+nVar*d+nchoosek(nVar, 2)*0.5*d*(d-1);
sMatrix = zeros(np,nVar);
coefVec = zeros(np,1);
if nargin < 4
   x0 = mean(vars.calBound,2);
else
   if ~vars.isFeasiblePoint(x0)
      x0 = mean(vars.calBound,2);
   end
end
if size(x0,1) ~= 1
   x0 = x0';
end

% 0-level HDMR
dm = x0;
coefVec(1) = f(dm);

% 1st-level HDMR
count = 1;
nSample = n1*d;
xS = vars.makeLHSsample(nSample);
for i = 1:nVar
   dm = repmat(x0, nSample, 1);
   dm(:,i) = xS(:,i);
   ytmp = f(dm);
   ytmp = ytmp - coefVec(1);
   xfit = expand1(xS(:,i) - x0(i));
   coefVec(count+1:count+d) = xfit\ytmp;
   sMatrix(count+1:count+d,i) = (1:d)';
   count = count+d;
end

% 2nd-level HDMR
nM = 0.5*d*(d-1);
nSample = n1*nM;
xS = vars.makeLHSsample(nSample);
tSM = zeros(nM,2);
cc = 1;
for i = 1:d-1 
   tSM(cc:cc+d-i-1,1) = i*ones(d-i,1);
   tSM(cc:cc+d-i-1,2) = (1:d-i)';
   cc = cc+d-i;
end
for i = 1:nVar-1
   for j = i+1:nVar
      dm = repmat(x0, nSample ,1);
      dm(:,[i, j]) = xS(:,[i, j]);
      Xi = expand1(xS(:,i) - x0(i));
      zi = (i-1)*d+1;
      fi = coefVec(zi+1:zi+d);
      yi = Xi*fi;
      Xj = expand1(xS(:,j) - x0(j));
      zj = (j-1)*d+1;
      fj = coefVec(zj+1:zj+d);
      yj = Xj*fj;
      ytmp = f(dm);
      ytmp = ytmp - coefVec(1)-yi-yj;
      xfit = expand2(xS(:,[i, j]) - repmat(x0([i,j]),nSample,1));
      coefVec(count+1:count+nM) = xfit\ytmp;
      sMatrix(count+1:count+nM,[i,j]) = tSM;
      count = count+nM;
   end
end
coefVec = reArrangeCoef(coefVec,x0);

pModel = B2BDC.B2Bmodels.PolyModel(sMatrix, coefVec, vars);

pModel.ErrorStats = estimateError(pModel);




   function xNew = expand1(xOld)
      nS = size(xOld,1);
      xNew = realpow(repmat(xOld,1,d),repmat(1:d,nS,1));
   end

   function xNew = expand2(xOld)
      xi = xOld(:,1);
      xj = xOld(:,2);
      xNew = zeros(size(xOld,1),0.5*d*(d-1));
      it = 0;
      for i1 = 1:d-1
         xNew(:,it+1:it+d-i1) = repmat(xi.^i1,1,d-i1).*...
            realpow(repmat(xj,1,d-i1),repmat(1:d-i1,length(xi),1));
         it = it+d-i1;
      end
   end

   function coefVec = reArrangeCoef(c,x0)
      if max(abs(x0)) > 1e-6
         coefVec = zeros(length(c),1);
         coefVec(1) = c(1);
         it = 1;
         for i1 = 1:nVar
            c1 = c(it+1:it+d);
            c2 = zeros(d,1);
            tmpx = x0(i1);
            for j1 = 1:d
               for j2 = j1:d
                  c2(j1) = c2(j1)+c1(j2)*nchoosek(j2,j1)*(-tmpx)^(j2-j1);
               end
               coefVec(1) = coefVec(1) + c1(j1)*(-tmpx)^j1;
            end
            coefVec(it+1:it+d) = c2;
            it = it+d;
         end
         n2 = 0.5*d*(d-1);
         for i1 = 1:nVar-1
            for i2 = i1+1:nVar
               xi = x0(i1);
               xj = x0(i2);
               idi = (i1-1)*d+1;
               idj = (i2-1)*d+1;
               c1 = c(it+1:it+n2);
               c2 = zeros(n2,1);
               for j1 = 1:n2
                  d1 = tSM(j1,1);
                  d2 = tSM(j1,2);
                  coefVec(1) = coefVec(1) + c1(j1)*(-xi)^d1*(-xj)^d2;
                  for z1 = 1:d1
                     coefVec(idi+z1) = coefVec(idi+z1) + c1(j1)*nchoosek(d1,z1)*(-xi)^(d1-z1)*(-xj)^d2;
                  end
                  for z2 = 1:d2
                     coefVec(idj+z2) = coefVec(idj+z2) + c1(j1)*nchoosek(d2,z2)*(-xj)^(d2-z2)*(-xi)^d1;
                  end
                  for k1 = 1:d1
                     for k2 = 1:d2
                        [~,tmpId] = intersect(tSM,[k1,k2],'rows');
                        c2(tmpId) = c2(tmpId) + c1(j1)*nchoosek(d1,k1)*(-xi)^(d1-k1)*nchoosek(d2,k2)*(-xj)^(d2-k2);
                     end
                  end
               end
               coefVec(it+1:it+n2) = c2;
               it = it+n2;
            end
         end
      else
         coefVec = c;
      end
   end

   function err = estimateError(pModel)
      nS = 2*np;
      xs = vars.makeLHSsample(nS);
      y1 = f(xs);
      y2 = pModel.eval(xs);
      dy = y1-y2;
      err.absMax = max(abs(dy));
      err.absAvg = mean(abs(dy));
      dy = abs(dy./y1);
      err.relMax = max(dy);
      err.relAvg = mean(dy);
   end


end