% The goal of this script is to replicate and verify the results of 
% "Comparison of statistic and deterministic frameworks of uncertainty 
% quantification" with the new B2B implementation.
% Recall, parameters in the paper are ordered with the old DCLABV2. 
% Thus, we need to order GRI variables as in DClabV2/CreateGRINew.

%% Load DCLabV2 GRI
D = createGRINew; %does not include the 37th QOI 'f5'
DCList = {D.FreeParameter(:).name}';

%% Load B2B GRI-Mech
[~,~,experimentData] = xlsread('GRIMech_expdata.xls');
s = load('GRIMech_modeldata.mat');
modelData = s.GRIMech_modeldata;
dsGRI = startDataset('GRI Mech 3.0');
dsGRI = addData(dsGRI,experimentData,modelData);
dsGRI.deleteUnit(37); %make the dataset consistent

B2BList = dsGRI.VarNames;
%% Check if the bounds for all QOIs are identical
% assuming the ordering is the same, for now. The naming convention between
% the two datasets is a bit different.
nQOI = dsGRI.Length;
violationIndex = [];
for kk = 1:nQOI
    dcQOI = D.ModelAndObservationPair(kk);
    b2bQOI = dsGRI.DatasetUnits.Values(kk);
    errLB = abs(dcQOI.constraintLowerLimit-b2bQOI.LowerBound);
    errUB = abs(dcQOI.constraintUpperLimit-b2bQOI.UpperBound);
    tol = 1e-2;
    if  errLB>tol || errUB>tol 
        violationIndex = [violationIndex; kk];
    end
end

% So, it does appear that there are differences between the QOI bounds (e.g. 'f6') in
% the two versions of the GRI-Mech dataset.
%% Find the order for the B2B variables
[intList, ind1, ind2] = intersect(DCList,B2BList,'stable');
assert(all(strcmp(DCList,B2BList(ind2)))==1)
%DCList == B2BList(ind2)

%% Test the prediction 
%Compare to Figure 1 in "Comparison of statistic and deterministic
%frameworks of uncertainty quantification". 

vars = dsGRI.Variables;
nVars = vars.Length;
Opt = B2BDC.Option({'Display',false});
for i1 = 1:numel(ind2)
    clc;
    disp(i1);
    xi = generateModel(vars.Values(ind2(i1)));
    dsGRI.setQOI2predict(xi);    
    dsGRI.predictQOI(Opt);
    bndsPosterior(i1) = dsGRI.QOIRange;
end

