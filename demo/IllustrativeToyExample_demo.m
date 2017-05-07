%% Illustrative Toy Example
% Simple three-component, 2-reaction linear model (3 state variables). 
% The experimental QOIs (eQOIs) are
%
% * peak-value of x2
% * time-of-peak-value of x2
% * time for x1 to decay to half initial value
%
% The prediction QOI (pQOI) is ratio of x3 to x1 at time of peak of x2.
% Model parameters are two rate coefficients


%% Two model parameters
% Create a |VariableList| object containing information about the model
% variables (uncertain/unknown rate constants).  More detailed information
% can be found in the VariableList_demo
varName = {'k1','k2'};
H = [1 4;0.5 2];
H = [1.3 3.6;0.6 1.8];
H = [1.3 3.6;0.6 1.8];
H = [[1.5 3];[1.5 3]];
varNominalVal = [2;1]
varNominalVal = [2.25;2.25]
varList = generateVar(varName,H,varNominalVal);

%% Create a random k data
nVar = size(H,1);
nSamples = 15*nVar*(nVar+1)/2;
kData = varList.makeLHSsample(nSamples);

%% Create Model 
% Set up 3-state ODE model for reaction system
odeModel = @(t,x,k1,k2) [-k1*x(1);k1*x(1)-k2*x(2);k2*x(2)];

%% Model simulation for Surrogate Model training
eSimData = zeros(nSamples,3);
pSimData = zeros(nSamples,1);
x0 = [1;0;0];
tFinal = 3;
for i=1:nSamples
   f = @(t,x) odeModel(t,x,kData(i,1),kData(i,2));
   [tSol,xSol] = ode45(f,[0 tFinal],x0);
   [x2peak,peakIndex] = max(xSol(:,2));
   [~,midIndex] = min(abs(xSol(:,1)-x0(1)/2));
   eSimData(i,:) = [x2peak tSol(peakIndex) xSol(midIndex,1)];
   pSimData(i,:) = xSol(peakIndex,3)/xSol(peakIndex,1);
end

%% Get quadratic fits to eQOIs
e1model = generateModelbyFit(kData,eSimData(:,1),varList,'q2norm');
e2model = generateModelbyFit(kData,eSimData(:,2),varList,'q2norm');
e3model = generateModelbyFit(kData,eSimData(:,3),varList,'q2norm');

%% Check accuracy of quadratic eQOI surrogate models
e1modelEvaluations = e1model.eval(kData,varList);
e2modelEvaluations = e2model.eval(kData,varList);
e3modelEvaluations = e3model.eval(kData,varList);
fitError1 = eSimData(:,1) - e1modelEvaluations;
fitError2 = eSimData(:,2) - e2modelEvaluations;
fitError3 = eSimData(:,3) - e3modelEvaluations;
meanError1 = mean(abs(fitError1)); stdError1 = std(fitError1);
meanError2 = mean(abs(fitError2)); stdError2 = std(fitError2);
meanError3 = mean(abs(fitError3)); stdError3 = std(fitError3);
disp(['e1QOI: Mean Abs Error: ' num2str(meanError1) ', STD Error: ' num2str(stdError1)]);
disp(['e2QOI: Mean Abs Error: ' num2str(meanError2) ', STD Error: ' num2str(stdError2)]);
disp(['e3QOI: Mean Abs Error: ' num2str(meanError3) ', STD Error: ' num2str(stdError3)]);

%% Get quadratic fits to pQOI
p1model = generateModelbyFit(kData,pSimData,varList);

%% Assess quality of quadratic surrogate model
p1modelEvaluations = p1model.eval(kData,varList);
fitError1 = pSimData - p1modelEvaluations;
meanError1 = mean(abs(fitError1)); stdError1 = std(fitError1);
disp(['p1QOI: Mean Abs Error: ' num2str(meanError1) ', STD Error: ' num2str(stdError1)]);

%% Rational quadratic fit to pQOI
p1modelrq = generateModelbyFit(kData,pSimData,varList,'rq');
p1modelrqEvaluations = p1modelrq.eval(kData,varList);
fitError1rq = pSimData - p1modelrqEvaluations;
meanError1rq = mean(abs(fitError1rq)); stdError1rq = std(fitError1rq);
disp(['p1QOIrq: Mean Abs Error: ' num2str(meanError1rq) ', STD Error: ' num2str(stdError1rq)]);

%% Synthetic Observed Values
% For the purposes of this example, we select a *true* value of the model
% variables, and compute the eQOIs at this true value.   Then,
% synthetically capture the experimental measurement by adding some
% plus/minus error to these true values to create the experimental bounds.
k1True = 2.1;
k2True = 2.4;
kTrue = [k1True k2True];
f = @(t,x) odeModel(t,x,k1True,k2True);
[tSol,xSol] = ode45(f,[0 tFinal],x0);
[x2peak,peakIndex] = max(xSol(:,2));
[~,midIndex] = min(abs(xSol(:,1)-x0(1)/2));
eTrue = [x2peak tSol(peakIndex) xSol(midIndex,1)];
pTrue = xSol(peakIndex,3)/xSol(peakIndex,1);

%%
% %% Check accuracy of quadratic eQOI surrogate models
% e1modelEvaluations = e1model.eval(kData,varList);
% e2modelEvaluations = e2model.eval(kData,varList);
% e3modelEvaluations = e3model.eval(kData,varList);
% fitError1 = eSimData(:,1) - e1modelEvaluations;
% fitError2 = eSimData(:,2) - e2modelEvaluations;
% fitError3 = eSimData(:,3) - e3modelEvaluations;
% meanError1 = mean(abs(fitError1)); stdError1 = std(fitError1);
% meanError2 = mean(abs(fitError2)); stdError2 = std(fitError2);
% meanError3 = mean(abs(fitError3)); stdError3 = std(fitError3);
% disp(['e1QOI: Mean Abs Error: ' num2str(meanError1) ', STD Error: ' num2str(stdError1)]);
% disp(['e2QOI: Mean Abs Error: ' num2str(meanError2) ', STD Error: ' num2str(stdError2)]);
% disp(['e3QOI: Mean Abs Error: ' num2str(meanError3) ', STD Error: ' num2str(stdError3)]);
% 
% 
% 
% %% Check accuracy of fit for Model #1
% % Plot the simulation and corresponding surrogate model predictions.  Note
% % that for this model, the surrogate model is quite accurate.
% e1modelEvaluations = e1model.eval(kData,varList);
% alle1 = [e1modelEvaluations; eSimData(:,1)];
% minmaxVal = [min(alle1) max(alle1)];
% plot(eSimData(:,1),e1modelEvaluations,'+',minmaxVal,minmaxVal);
% xlabel('Simulation Data: #1')
% ylabel('(quadratic) Surrogate Model #1 prediction');
% title('Comparing model to surrogate model')
% fitError1 = eSimData(:,1) - e1modelEvaluations;
% meanError1 = mean(abs(fitError1));
% stdError1 = std(fitError1);
% disp(['Mean Abs Error: ' num2str(meanError1) ', STD Error: ' num2str(stdError1)]);
% 
% 
% 
% %% Check accuracy of fit for Model #2
% % Plot the simulation and corresponding surrogate model predictions.  
% e2modelEvaluations = e2model.eval(kData,varList);
% alle2 = [e2modelEvaluations; eSimData(:,2)];
% minmaxVal = [min(alle2) max(alle2)];
% plot(eSimData(:,2),e2modelEvaluations,'+',minmaxVal,minmaxVal);
% xlabel('Simulation Data: #2')
% ylabel('(quadratic) Surrogate Model #2 prediction');
% title('Comparing model to surrogate model')
% fitError2 = eSimData(:,2) - e2modelEvaluations;
% meanError2 = mean(abs(fitError2));
% stdError2 = std(fitError2);
% disp(['Mean Abs Error: ' num2str(meanError2) ', STD Error: ' num2str(stdError2)]);
% 
% %% Rational-quadratic for Model #2
% e2modelrq = generateModelbyFit(kData,eSimData(:,2),varList,false,false,'rq');
% e2modelrqEvaluations = e2modelrq.eval(kData,varList);
% alle2 = [e2modelrqEvaluations; eSimData(:,2)];
% minmaxVal = [min(alle2) max(alle2)];
% plot(eSimData(:,2),e2modelEvaluations,'+',...
%    eSimData(:,2),y2modelrqEvaluations,'go',...
%    minmaxVal,minmaxVal);
% xlabel('Simulation Data')
% ylabel('(quadratic and rational-quadratic) Surrogate Model prediction');
% title('Comparing model to surrogate model')
% fitError2rq = eSimData(:,2) - y2modelrqEvaluations;
% meanError2rq = mean(abs(fitError2rq));
% stdError2rq = std(fitError2rq);
% disp(['Mean Abs Error: ' num2str(meanError2rq) ', STD Error: ' num2str(stdError2rq)]);
% 
% %% Get quadratic fit to pQOI
% p1model = generateModelbyFit(kData,pSimData,varList);
% 
% %% Assess quality of quadratic surrogate model
% p1modelEvaluations = p1model.eval(kData,varList);
% allp1 = [p1modelEvaluations; pSimData];
% minmaxVal = [min(allp1) max(allp1)];
% plot(pSimData,p1modelEvaluations,'+',minmaxVal,minmaxVal);
% xlabel('Simulation Data: Prediction Model')
% ylabel('(quadratic) Surrogate Model');
% title('Comparing model to surrogate model')
% fitError1 = pSimData - p1modelEvaluations;
% meanError1 = mean(abs(fitError1));
% stdError1 = std(fitError1);
% disp(['Mean Abs Error: ' num2str(meanError1) ', STD Error: ' num2str(stdError1)]);
% 
% 
%% Experimental Observations regarding eQOI1 and eQOI2
% These are bounds, asserted by experimentalist, as to the true values
% (from experimentation) of the two experimental QOIs.
eBound1 = [0.55 0.83];
eBound2 = [0.45 0.55];
eBound3 = [0.45 0.55];
eBound1 = eTrue(1)*[0.9 1.1];
eBound2 = eTrue(2)*[0.9 1.1];
eBound3 = eTrue(3)*[0.9 1.1];
% 
% %% Create DataSet Units for all eQOIs
% % Each DataSet Unit is associated with an experimental QOI, and consists of
% % a (surrogate) model and associated bounds on the eQOI's true value. 
dsUnitName1 = 'Peak X2';
dsUnit1 = generateDSunit(dsUnitName1,e1model,eBound1);
dsUnitName2 = 'Time of Peak X2';
dsUnit2 = generateDSunit(dsUnitName2,e2model,eBound2);
dsUnitName3 = 'Time of Half X1';
dsUnit3 = generateDSunit(dsUnitName3,e3model,eBound3);
% 
%%  Create a DataSet object, and insert all DataSet Units
dsName = 'Toy SIAM Example';
DS = generateDataset(dsName);
DS.addDSunit(dsUnit1);
DS.addDSunit(dsUnit2);
DS.addDSunit(dsUnit3);

%% Consistency Test: check consistency of DataSet
DS.isConsistent


%% Prediction
[pRange, xFeasInner, Sens] = DS.predictQOI(p1model);
% 
% 
% 
% 
% 
