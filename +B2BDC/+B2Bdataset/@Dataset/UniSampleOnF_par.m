function [xVal,eff] = UniSampleOnF_par(obj, N, xS, Dinfo, opt)
% XVAL = UNISAMPLEONF_PAR(OBJ,N,XS,DINFO,OPT) returns a N-by-nvariable sample matrix.
% Each row of the matrix corresponds to a sample point uniformly drawn from
% the feasible set of the obj dataset. 

%  obj - A B2BDC.B2Bdataset.Dataset object
%  N - Number of samples returned
%  xS - nSample-by-nVar RW sample matrix for principal direction analysis
%       or uncertainty bound calculation. If empty, the code will generate
%       RW samples automatically if needed
% Dinfo - A structure variable with fields D and uq. For RSH method the D
%         specifies the direction of the box by the column vectors and uq specifies
%         the span along each direction. For RSP method, D specifies the
%         directions for the polytope and uq specifies the bounds along
%         each direction. If it is empty, RSH will use sample estimated PC
%         directions and RSP will use extra directions specified by
%         opt.SampleOption.ExtraCut. The direction matrix is specified by
%         the column of D.
% opt - A B2BDC.Option.Option object

%  Created: May 2, 2016

if nargin < 5
   opt = generateOpt('Display',false);
end
if nargin < 3
   xS = [];
end
if nargin < 4
   Dinfo = [];
end
sopt = opt.SampleOption;
if ~obj.isConsistent(opt)
   disp('The dataset is not consistent')
   xVal = [];
else
   smethod = sopt.SampleMethod;
   switch(smethod)
      case 'RSH'
         [xVal,eff] = obj.RSH_par(N,xS,Dinfo,opt);
      case 'RSP'
         [xVal,eff] = obj.RSP_par(N,xS,Dinfo,opt);
   end
end