function polyBasis = generatePolyBasis(pBasis)
% Create a B2BDC.B2Bvariables.PolyBasis object.
% The input argument is:
%   pBasis - A nMonomial-by-nVariable matrix specifies the basis

%  Created: Nov 15, 2015     Wenyu LiB
if ismatrix(pBasis)
   polyBasis = B2BDC.B2Bvariables.PolyBasis(pBasis);
else
   error('The input basis should be a nBasis-by-nVar matrix')
end
