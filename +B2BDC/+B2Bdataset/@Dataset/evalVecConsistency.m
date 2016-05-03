function [newObj, vecCM] = evalVecConsistency(obj,W)
%%EVALVECCONSISTENCY
%prototype vector consistency measure for use with GRImech
% Arun Hegde 10/5/2015

%% Step 1: Store relevant information from dataset obj
numVar =  numel(obj.Variables.Values); %number of variables
numExp = numel(obj.DatasetUnits.Values); %number of experiments
Bnds = [[obj.DatasetUnits.Values.LowerBound]' [obj.DatasetUnits.Values.UpperBound]']; %[lb ub measurement];
if nargin == 1
    W = diag(ones(numExp));
end
%% Step 2: Build relevant objective and constraint matrices
% Terminology comes from formulations presented in powerpoint and pdf document.
Cobj = zeros(1+numVar+numExp+numExp,1+numVar+numExp+numExp);
Qu = cell(numExp,1);
Ql = cell(numExp,1);
Pu = cell(numExp,1);
Pl = cell(numExp,1);
H = cell(numVar,1);

I = eye(numExp,numExp); 
for jj = 1:numExp
    M = expandModels(obj, jj);   %expands the jjth model from active variables to full variables.
    A = M(2:numVar+1, 2:numVar+1);
    c = M(1,1);
    b = M(2:numVar+1,1);
    u = Bnds(jj,2);
    l = Bnds(jj,1);
    
    Qu{jj} = [c-u b' zeros(1,numExp) -0.5*I(jj,:) ; ...
        b A zeros(numVar,numExp) zeros(numVar,numExp) ; ...
        zeros(numExp,1) zeros(numExp,numVar) zeros(numExp,numExp) zeros(numExp, numExp); ...
        -0.5*I(:,jj) zeros(numExp,numVar) zeros(numExp,numExp) zeros(numExp, numExp)];
    
    Ql{jj} = [-c+l -b'  -0.5*I(jj,:) zeros(1,numExp); ...
        -b -A zeros(numVar,numExp) zeros(numVar,numExp) ; ...
        -0.5*I(:,jj) zeros(numExp,numVar) zeros(numExp,numExp) zeros(numExp, numExp); ...
        zeros(numExp,1) zeros(numExp,numVar) zeros(numExp,numExp) zeros(numExp, numExp)];
    
    Pu{jj} = [0 zeros(1,numVar+numExp) -0.5*I(jj,:); ...
        zeros(numVar,1) zeros(numVar,numVar+numExp+numExp) ; ...
        zeros(numExp,1+numVar+numExp+numExp); ...
        -0.5*I(:,jj) zeros(numExp,numVar+numExp+numExp)];
    
    Pl{jj} = [0 zeros(1,numVar) -0.5*I(jj,:) zeros(1,numExp); ...
        zeros(numVar,1) zeros(numVar,numVar+numExp+numExp) ; ...
        -0.5*I(:,jj) zeros(numExp,numVar+numExp+numExp); ...
        zeros(numExp,1+numVar+numExp+numExp)];
end

for jj = 1:numVar  %bounds on X
    Ix = zeros(numVar,numVar);
    Ix(jj,jj) = 1;
    H{jj} = [-1 zeros(1,numVar+numExp+numExp); ...
        zeros(numVar,1) Ix zeros(numVar,numExp+numExp); ...
        zeros(numExp+numExp,1+numVar+numExp+numExp)];
end

constrSet = [Qu; Ql; Pu; Pl; H];

Cobj = [0 zeros(1,numVar) 0.5*ones(1,numExp) 0.5*ones(1,numExp); ...
    zeros(numVar,1) zeros(numVar,numVar) zeros(numVar,numExp) zeros(numVar,numExp); ...
    0.5*ones(numExp,1) zeros(numExp,numVar) zeros(numExp,numExp) zeros(numExp, numExp); ...
    0.5*ones(numExp,1) zeros(numExp,numVar) zeros(numExp,numExp) zeros(numExp, numExp)];

Wmat = diag(W);
Wfull = blkdiag(eye(numVar+1),Wmat,Wmat);
Cobj = Wfull'*Cobj*Wfull;

%% Step 3: Sedumi implementation of the optimization in primal form & formatting sedumi result
% Details regarding derivation can be found in "A note on the quadratic vector consistency
% measure: converting primal form to sedumi" draft pdf document.
% x = [slackvar; vec(Z)];
numConstr = numel(constrSet);
p = numel(Cobj);
c = [zeros(numConstr,1); vec(Cobj')];
A1 = zeros(numConstr, p);
for kk = 1:numConstr
    A1(kk,:) = vec(constrSet{kk}')';
end
A = [eye(numConstr,numConstr) A1; ...
    zeros(1,numConstr) 1 zeros(1,p-1)];
b = [zeros(numConstr,1); 1];

K.l = [numConstr];
K.s = [sqrt(p)];

[x,y, info] = sedumi(A,b,c,K);

Z = mat(x(numConstr+1:end));
z = Z(2:end,1);
E = Z(2:end,2:end);


%% Step 5: Finding an upper bound via fmincon
%
% cvx_begin SDP
%
% cvx_solver sedumi
% %variable E(numVar+numExp+numExt,numVar+numExp+numExp) symmetric
% %variable z(numVar+numExp+numExp,1) % z = [x dl du]';
% variable EE(size(Cobj)-1) symmetric
% variable k(size(Cobj,1)-1,1) % z = [x dl du]';
%
% dual variable y{numConstr};
% K = [1 k'; k EE];
%
% minimize trace(Cobj*Z)
% subject to
% K >= 0;
% for jj = 1:numConstr
%     trace(constrSet{jj}*K) <= 0 : y{jj};
% end
%
% cvx_end
CovZ = E-z*z';
func = @(x) objFn(x,Cobj);
ineqC = @(x) constrFn(x,constrSet);
options = optimset('Display', 'on','GradObj','on','GradConstr','on','Hessian','user-supplied', ...
    'HessFcn',@(x,lambda)quadhess(x,lambda,Cobj,constrSet)) ;
[zFeas objFeas] = fmincon(func,z,[],[],[],[],-inf*ones(size(Cobj,1)-1,1),inf*ones(size(Cobj,1)-1,1), ineqC, options);

vecCM.name = 'Vector consistency measure';
vecCM.lb = c'*x;
vecCM.ub = objFeas;
vecCM.x = zFeas(1:numVar);
vecCM.dl = zFeas(numVar+1:numVar+numExp);
vecCM.du = zFeas(numVar+numExp+1:numVar+numExp+numExp);
vecCM.Id = unique([find(vecCM.dl>1e-4);find(vecCM.du>1e-4)]);

% figure;
% subplot(1,2,1);
% bar(vecCM.dl);
% title('Lower bound relaxation');
% xlabel('Constraint Index');
% ylabel('Magnitude of relaxation');
% subplot(1,2,2);
% bar(vecCM.du);
% title('Upper bound relaxation');
% xlabel('Constraint Index');
% ylabel('Magnitude of relaxation');

%% new dataset object with the adjusted bounds. Finish by verifying consistentcy via the old scalar metric. 
newObj = clone(obj);
newLB = Bnds(:,1) - vecCM.dl;
newUB = Bnds(:,2) + vecCM.du;
newObs = [newObj.DatasetUnits.Values.ObservedValue]';
newObj.changeBounds([newLB, newUB,newObs],[1:numExp]');
newObj.FeasiblePoint = zFeas(1:numVar);
% newObj.isConsistent
end

%% Extra helper functions

function [Q] = expandModels(dsObj, ind)
%For use with GRI Mech 3.0

varList = {dsObj.Variables.Values.Name}';
activeVar = {dsObj.DatasetUnits.Values(ind).VariableList.Values.Name}';

nAct = numel(activeVar);
nVar = numel(varList);
index = zeros(nAct,1);
Q = zeros(nVar+1,nVar+1);

for kk = 1:nAct
    index(kk) = find(strcmp(activeVar{kk}, varList));
end

index1 = [1; index+1];

Q(index1, index1) = dsObj.DatasetUnits.Values(ind).SurrogateModel.CoefMatrix;

end


function [obj, Grad] = objFn(x, Cobj)

obj = [1;x]'*Cobj*[1;x];
nX = size(Cobj,1) -1;
A = Cobj(2:nX+1, 2:nX+1);
b = 2*Cobj(2:nX+1,1);

if nargout > 1
    Grad = (A+A')*x+b;
end
end

function [ineqConstr, eqConstr, gradIneqConstr, gradEqConstr] = constrFn(x, Constraints)

ineqConstr = zeros(numel(Constraints),1);
nX = size(Constraints{1},1)-1;
if nargout >2
    gradIneqConstr = zeros(nX,numel(Constraints));
end
for jj = 1:numel(Constraints)
    ineqConstr(jj,1) = [1;x]'*Constraints{jj}*[1;x];  %related to the original problem, NOT the relaxed SDP
    
    if nargout > 2
        A = Constraints{jj}(2:nX+1, 2:nX+1);
        b = 2*Constraints{jj}(2:nX+1,1);
        gradIneqConstr(:,jj) = (A+A')*x +b;
    end
end

eqConstr = [];
gradEqConstr = [];
end

function hess = quadhess(x,lambda,Cobj,Constraints)
nX = size(Cobj,1) -1;
Ao = Cobj(2:nX+1, 2:nX+1);
hess = Ao+Ao';
n = numel(Constraints);
for ii = 1:n
    Ac = Constraints{ii}(2:nX+1, 2:nX+1);
    hess = hess + lambda.ineqnonlin(ii)*(Ac+Ac');
end
end
