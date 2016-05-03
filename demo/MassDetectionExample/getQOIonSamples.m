function ySamples = getQOIonSamples(xSamples, wVec, nMass, fixedMassVec, Kval, Cval)

nMassFixed = numel(fixedMassVec);
szX = size(xSamples);
szw = size(wVec);
if numel(szX)==2 && isvector(wVec)
   nS = szX(1);
   nX = szX(2);
   N = nMassFixed + nX;
   nY = numel(wVec);
   ySamples = zeros(nS,nY);
   for i=1:nS
      massValue = [xSamples(i,:)'; fixedMassVec];
      [A,B,C] = msd1d(massValue,Kval,Cval);
      z = zeros(nY,1);
      for p = 1:nY
      H = (C*(wVec(p)*1i*eye(2*N) - A)\B);
      z(p) = H(N+1,1);
      end
      ySamples(i,:) = abs(z(:))';
   end
else
    error('xSamples should be 2-d, nSamples-by-nX; wVec should be a vector');
end
