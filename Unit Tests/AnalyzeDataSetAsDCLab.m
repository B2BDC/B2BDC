
%% Ryan/Trent DCLab DataSet
clear
load Friday18December

%%
ntargets = numel(MatsVars);
AllVars = cell(0,1);
trgVars = cell(ntargets,1);
trgNames = cell(ntargets,1);
for i=1:ntargets
   trgVars{i} = MatsVars(i).Vars;
   AllVars = union(AllVars,trgVars{i});
end
nX = numel(AllVars);

Vidx = cell(ntargets,1);
for i=1:ntargets
   [~,b,c] = intersect(AllVars,trgVars{i});
   % AllVars(b) equals trgVars{i}(c), and "c" reshuffles trgVars{i}.
   % want AllVars(Vidx{i}) to equal trgVars{i} (ie., unschuffled)
   [~,idx] = sort(c);
   Vidx{i} = b(idx);
   assert(isequal(AllVars(Vidx{i}),trgVars{i}),...
     ['Target ' int2str(i) ' variables and AllVars are not consistent']);
end

Mfull = zeros(1+nX,1+nX,ntargets);
Mcell = cell(ntargets,1);  % only contains active entries
UB = zeros(ntargets,1);
LB = zeros(ntargets,1);
for i=1:ntargets
   Mcell{i} = MatsVars(i).Matrix;
   idx = [1;1+Vidx{i}];
   Mfull(idx,idx,i) = MatsVars(i).Matrix;
   UB(i) = MatsVars(i).YBounds(2);
   LB(i) = MatsVars(i).YBounds(1);
end

%% Create DataSet
nparms = length(AllVars);
model_domain = struct('name',cell(nparms,1),'range',cell(nparms,1));
FP = cell(nparms,1);
for i=1:nparms
   model_domain(i).name = AllVars{i};
   model_domain(i).range = [-10 10];
   FP{i} = FreeParameter(AllVars{i},0,[-1 1]);
end
RM = cell(ntargets,1);
OP = cell(ntargets,1);
RMOP = cell(ntargets,1);
for i=1:ntargets
   RM{i} = ResponseModel(Mcell{i},model_domain(Vidx{i}));
   Ysyn = 0.5*(UB(i) + LB(i));
   R = 0.5*(UB(i)-LB(i));
   OP{i} = ResponseObservation(Ysyn,R);
   RMOP{i} = ModelAndObservationPair(OP{i},RM{i},['Y' int2str(i)]);
end
MOPArray(1) = RMOP{1};
for i=2:ntargets
   MOPArray(i) = RMOP{i};
end
FPArray(1) = FP{1};
for i=2:nparms
   FPArray(i) = FP{i};
end
D = DCDataset(MOPArray,FPArray);

%% Consistency Test of DataSet
% The ConsistencyTest maximizes the amount each dataset constraint can be
% tightened, under the restriction that the DataSet remain feasible.  If
% the DataSet is originally infeasible, then the consistency measure will
% be negative (indicating that constraints must be relaxed to obtain
% feasibliity).  The consistency test is a *maximization*, so the upper
% bound is an outer bound, and the lower bound is an inner bound.
opts = DCOptions;
opts.nRestart = 2;
opts.maxBranchBoundIter = 0;
C = ConsistencyTest(D,opts);

%% Posterior predictions: Parameters
nX = nparms;
RPs = struct('UBo',[],'LBo',[],'UBi',[],'LBi',[],'UBx',[],'LBx',[]);
RPs = repmat(RPs,[nX 1]);
opts = DCOptions;
opts.maxBranchBoundIter = 1;
for i=1:nX
   xRM = ResponseModel([0 .5;.5 0],model_domain(i));
   tic
   RP = ResponsePrediction(xRM,D,opts);
   [i toc]
   RPs(i).UBo = RP.UBo;
   RPs(i).LBo = RP.LBo;
   RPs(i).UBi = RP.UBi;
   RPs(i).LBi = RP.LBi;
   RPs(i).UBx = RP.UBx;
   RPs(i).LBx = RP.LBx;
end



