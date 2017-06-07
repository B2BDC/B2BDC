%% B2BDC example: how to use various sampling routines
% In this demo, examples of how to use all the sampling routines (random
% walk, rejection with box, rejection with polytope) are given. Examples of 
% optional utilization of approximation strategy as well as dimension reduction 
% through principal component analysis are also given.

%% create a feasible dataset
ds = generateDataset('Example dataset for sampling');
nVar = 5;
nQ = 10;
tVar = generateVar([],repmat([-1,1],nVar,1));
for i = 1:nQ
   CoefMatrix = rand(nVar+1);
   CoefMatrix(1,2:end) = 0.5*CoefMatrix(1,2:end); % scaling linear part
   CoefMatrix(2:end,1) = 0.5*CoefMatrix(2:end,1);
   CoefMatrix(2:end,2:end) = 0.1*CoefMatrix(2:end,2:end); % scaling quadratic part
   CoefMatrix = 0.5*(CoefMatrix+CoefMatrix');
   lb = CoefMatrix(1,1) - rand(1);
   ub = CoefMatrix(1,1) + rand(1);
   tmp_Model = generateModel(CoefMatrix,tVar);
   tmp_DSunit = generateDSunit(['Unit ' num2str(i)], tmp_Model,[lb,ub]);
   ds.addDSunit(tmp_DSunit);
end

%% generate B2BDC.Option object
opt = generateOpt('Display',false);

%% Check the feasibility of the dataset (it should be consistent)
ds.isConsistent(opt)

%% All sampling-related options are in opt.SampleOption
opt.SampleOption

%% generate samples in the feasible set by random walk
% The only controlling parameter in this method is the filed 'StepInterval'
% which specifies the algorithm to take only the k th point in the chain.
% k=1 corresponds to the situation that the whole chain is taken as the
% final sample.
nS = 1e3;
x_start = [];
x_rw = ds.collectSamples(nS,[],opt);

% You can also change the step interval or the starting point by setting
% them before passing into the function
opt.SampleOption.StepInterval = 20;
x_start = ds.FeasiblePoint;  % this point has to be a feasible point
x_rw = ds.collectSamples(nS,[],opt);

%% generate uniform samples of the feasible set with bounding box
% set sampling method to RSH (rejection sampling with hyperrectangle)
opt.SampleOption.SampleMethod = 'RSH';
% set the batch size of each candidate generation
opt.SampleOption.BatchMaxSample = 1e5;
% first use the original box H
Dinfo.D = eye(nVar);
Dinfo.uq = ds.Variables.calBound;
% The detailed input explaination can be found in the function script
nS = 1e3;
xF = []; % no need for rw samples in this case
[x_rsh,eff_rsh] = ds.UniSampleOnF(nS,xF,Dinfo,opt);
eff_rsh % the calculated sampling efficiency

% now use the box defined by principal directions, the uncertainy bound is
% using the outer bound estimation (no approximation)
Dinfo = [];
opt.SampleOption.UncertaintyEstimation = 'Outer';
[x_rsh,eff_rsh] = ds.UniSampleOnF(nS,xF,Dinfo,opt);
eff_rsh

%% generate uniform samples of the feasible set with polytope
% specifying the direction and uncertainty manually
nL = 50; % number of directions specifying the polytope
D = randn(nVar,nL);  % generate random directions
D = normc(D);
uq = zeros(nL,2);
opt.Prediction = 'outer'; % using the outer bound estimation
% calculate the uncertainy bounds
for i = 1:nL
   tModel = generateModel([0;D(:,i)],tVar); 
   p = ds.predictQOI(tModel,opt);
   uq(i,:) = [p.min, p.max]; % using the outer bound estimation
end
% create the Dinfo input
Dinfo.D = D;
Dinfo.uq = uq;
% set sampling method to RSP
opt.SampleOption.SampleMethod = 'RSP';
% set xF empty since there is no need in this case
xF = [];
nS = 1e3;
[x_rsp, eff_rsp] = ds.UniSampleOnF(nS,xF,Dinfo,opt);
% the sampling efficiency
eff_rsp

%% approximation strategy
% approximation strategy is included by using a polytope containing
% trunaced feasible set

% #1 using inner bound instead of outer bound to estimate the polytope
%    bounds
opt.SampleOption.UncertaintyEstimation = 'Inner';
% instead of create D mannually, using the sample option
Dinfo = [];
opt.SampleOption.ExtraCut = 50;
% set xF empty since there is no need in this case
xF = [];
nS = 1e3;
% generate samples
[x_rsp, eff_rsp] = ds.UniSampleOnF(nS,xF,Dinfo,opt);
% the sampling efficiency
eff_rsp

% #2 using sample bound instead of outer bound to estimate the polytope
%    bounds
opt.SampleOption.UncertaintyEstimation = 'Sample';
% instead of create D mannually, using the sample option
opt.SampleOption.ExtraCut = 50;
Dinfo = [];
% generate xF for the sampling estimation invloved in this case
xStart = [];
xF = ds.collectSamples(1e4,xStart,opt);
% generate samples
[x_rsp, eff_rsp] = ds.UniSampleOnF(nS,xF,Dinfo,opt);
% the sampling efficiency
eff_rsp

%% dimension reduction with PCA
% trucate 2 principal directions
opt.SampleOption.TruncatedPC = 2;
% using outer bound estimation
opt.SampleOption.UncertaintyEstimation = 'Outer';
opt.SampleOption.ExtraCut = 50;
% generate xF for the PCA invloved in this case
xStart = [];
xF = ds.collectSamples(1e4,xStart,opt);
% generate samples
nS = 1e3;
Dinfo = [];
[x_rsp, eff_rsp] = ds.UniSampleOnF(nS,xF,Dinfo,opt);
% the sampling efficiency
eff_rsp

%% parallel version of the sampling code
% simply replace the call ds.UniSampleOnF by ds.UniSampleOnF_par for all
% the cases shown above

%% sampling a polytope
% assume you start with a prior domain H
nVar = 10;
H = [-rand(nVar,1) rand(nVar,1)];
varList = generateVar([],H);
% assume you find other linear constraints such that lb < Ax <ub
nExtra = 20;
A = randn(nExtra,nVar);
lb = -rand(nExtra,1);
ub = rand(nExtra,1);
% add the linear constraints to varList so it represents the polytope
varList_polytope = varList.addLinearConstraint(A,lb,ub);
% set start point to be random
xStart = [];
% set step size for a better mixing
opt.SampleOption.StepInterval = 20;
% collect random walk samples in the polytope (which is uniformly
% distributed)
nS = 1e3;
x_rw = varList_polytope.collectSamples(nS,xStart,opt.SampleOption);