%% B2BDC example: Using vector consistency to analyze an inconsistent dataset
% In this example, we demonstrate a vector consistency analysis of the 
% GRI-Mech 3.0 dataset, which was shown to be inconsistent in |GRIMech_demo2.m|.    

%% Loading the GRI-Mech 3.0 dataset 
% The GRI-Mech 3.0 dataset |dsGRI| is constructed from an Excel 
% spreadsheet and a MATLAB .mat file. For a detailed walkthrough of these
% steps, please refer to |GRIMech_demo1.m|.

[~,~,experimentData] = xlsread('GRIMech_expdata.xls');
s1 = load('GRIMech_modeldata.mat');
modelData = s1.GRIMech_modeldata;
dsGRI = generateDataset('GRI Mech 3.0');
dsGRI = addData(dsGRI,experimentData,modelData);

%% 
% In this example, the vector consistency measure (VCM) is applied to |dsGRI| with
% 'unit' weights on the QOI bounds and 'null' weights on the variable bounds. 
% The output structure, |vcReport1| is displayed below.

vcReport1 = dsGRI.vectorConsistency('unit','null');
vcReport1

%%
% |vcReport1.Objective| provides bounds on the objective of the VCM. 
% The lower bound is calculated by a semidefinite program relaxation 
% and the upper bound (a local solution) comes from |fmincon|.
% |vcReport1.FeasiblePoint| stores a point which is feasible
% if the relaxations in |vcReport1.Relaxations| are implemented. 
% |vcReport1.Weights| stores the user-specified weights
% related to the relaxations. When implementing a relaxation, the
% values in |vcReport1.Relaxations| must be scaled by the corresponding weights. 

vcReport1.Objective
vcReport1.Relaxations
vcReport1.Weights

%% 
% Examining |vcReport1| reveals that consistency can be reached by relaxing  
% the upper bound of QOI #36 and the lower bound of QOI #37. 

relLBIdx = find(vcReport1.Relaxations.yLowerBound > 1e-5);
relUBIdx = find(vcReport1.Relaxations.yUpperBound > 1e-5);
relIdx = union(relLBIdx,relUBIdx);
numel(relIdx);

relUBIdx
relLBIdx



