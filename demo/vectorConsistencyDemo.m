%% B2BDC example: Using vector consistency to analyze inconsistent datasets
% In this example, we demonstrate a vector consistency analysis of the 
% GRI-Mech 3.0 dataset and the DLR-SynG dataset. 

%% Loading the GRI-Mech 3.0 dataset 
% For a detailed walkthrough of these steps, please refer to |GRIMech_demo1.m|.

[~,~,experimentData] = xlsread('GRIMech_expdata.xls');
s1 = load('GRIMech_modeldata.mat');
modelData = s1.GRIMech_modeldata;
dsGRI = generateDataset('GRI Mech 3.0');
dsGRI = addData(dsGRI,experimentData,modelData);


%%
% In |GRIMech_demo2.m|, the notion of *dataset consistency* is defined and 
% |dsGRI| is shown to be inconsistent by computing the consistency measure.
% The vector consistency measure (VCM) is an additional tool for analyzing
% inconsistent datasets that allows independent relaxations to each
% model-data and prior constraint, i.e.
%
% $$L_e - W_e^{L} \Delta_{e}^{L} \leq M_e(x) \leq U_e + W_e^{U}
%  \Delta_{e}^{U}$$
%
% for e indexing the model-data constraints (QOIs) and
%
% $$l_i - w_i^{l} \delta_{i}^{l} \leq x_i \leq u_i + w_i^{u}
%  \delta_{i}^{u}$$
% 
% where i indexes the variables. The VCM minimizes the 1-norm sum of the relaxations, i.e.
%
% $$\min \|\Delta^L\|_1 + \|\Delta^U\|_1 + \|\delta^l\|_1 + \|\delta^u\|_1$$
% 
% subject to the above constraints. The 1-norm is a well-known heuristic 
% for sparsity.  The coefficients to the relaxations (W,w) act as weights, 
% where smaller values protect against relaxation and larger values promote 
% relaxation to the corresponding constraint. For instance, 
%
% $$W_3^{L} = 0$$ 
%
% prevents relaxation to the lower bound of the 3rd QOI, thus reinforcing 
% the original experimental observations. 

%%
% The method |vectorConsistency| computes the VCM for a given dataset. The
% help for the |vectorConsistency.m| is shown below.

help vectorConsistency

%% 
% In this section, the vector consistency measure (VCM) is applied to |dsGRI| with
% 'unit' weights on the QOI bounds and 'null' weights on the variable bounds. 
% The output structure, |vcReport1| is displayed below.

vcReport1 = dsGRI.vectorConsistency('unit','null', 5);
vcReport1

%%
% |vcReport1.Objective| provides bounds on the objective of the VCM, thus 
% bracketing the true value. The lower bound is calculated by a semidefinite program (SDP)  
% and the upper bound (a local solution) comes from |fmincon|. The SDP
% result tells us that no set of relaxations that sums to less than 0.0174
% can lead to consistency. The local solution provides a collection of
% relaxations that are feasible, i.e. that if implemented with the specified weights 
% produce a consistent dataset.

vcReport1.Objective

%%
% |vcReport1.FeasiblePoint| stores a point which is feasible
% if the relaxations in |vcReport1.Relaxations| are implemented with the 
% specified weights. The feasible point and the relaxations come from the 
% local solution. 

vcReport1.FeasiblePoint;
vcReport1.Relaxations

%%
% |vcReport1.Weights| stores the user-specified weights
% related to the relaxations.

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

%%
% Utilizing the |perc| weighting scheme, detailed in the VCM help, reveals that
% the GRI-Mech dataset can be made consistent by adjust just a single
% bound. In fact, the dataset can be made consistent by decreasing the
% lower bound of QOI 37 by just under 1.6%. Hence, the inconsistency of
% |dsGRI| can be explained by a single QOI.

vcReport2 = dsGRI.vectorConsistency('perc','null',1);
relLBIdx = find(vcReport2.Relaxations.yLowerBound > 1e-5);
relUBIdx = find(vcReport2.Relaxations.yUpperBound > 1e-5);
vcReport2.Relaxations.yLowerBound(relLBIdx)

%%
% In this section, we switch gears and investigate the DLR-SynG dataset, 
% labeled as |dsDLR| dataset. This dataset is comprised of 159 QOIs and
% 55 uncertain parameters. As will be illustrated, the degree of 
% inconsistency in |dsDLR| is more pronounced as compared to |dsGRI|. 

load dsDLR_SynG.mat
dsDLR.DatasetUnits
dsDLR.Variables

%%
% First, note that the dataset is inconsistent. Fitting error is included
% in the analysis.

Opt = generateOpt('Display',false,'AddFitError',true);
dsDLR.isConsistent(Opt)
dsDLR

%%
% The VCM is now applied to |dsDLR| with 'unit' weights on the QOI bounds 
% and 'null' weights on the variable bounds. The output structure is 
% displayed below, as is the number of constraints relaxed in the local solution. 

vcReport3 = dsDLR.vectorConsistency('unit','null',5, Opt)
vcReport3.Objective
relLBIdx = find(vcReport3.Relaxations.yLowerBound > 1e-5);
relUBIdx = find(vcReport3.Relaxations.yUpperBound > 1e-5);
relIdx = union(relLBIdx,relUBIdx);
numel(relIdx)

%%
% Different weighting schemes can result in fewer constraint relaxations. 
% Here, we use the 'uwidth' weighting configuration. 

vcReport4 = dsDLR.vectorConsistency('uwidth','null',5, Opt)
vcReport4.Objective
relLBIdx = find(vcReport4.Relaxations.yLowerBound > 1e-5);
relUBIdx = find(vcReport4.Relaxations.yUpperBound > 1e-5);
relIdx = union(relLBIdx,relUBIdx);
numel(relIdx)

%% Attribution
% UC Berkeley, Spring 2017, B2BDC team