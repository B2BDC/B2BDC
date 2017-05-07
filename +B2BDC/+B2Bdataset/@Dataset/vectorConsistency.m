function [vcReport, EXITFLAG] = vectorConsistency(dsObj,wY,wX, nInit, opt)
% vectorConsistency   Evaluates the vector consistency of a given dataset
% with given weights
%    V = obj.vectorConsistency(wY,wX)  returns a structure containing the
%    results of the vector consistency analysis, including the necessary
%    relaxations to the QOI bounds and the variable bounds to render the
%    dataset consistent.
%
% The inputs are:
% --------------------------------------------------------------------------
%     dsObj - An inconsistent dataset
%     wY    - A obj.Length-by-2 matrix of weights, where the first column
%            acts on the relaxations to the corresponding QOI lower bounds
%            and the second column acts on the QOI upper bounds. A weight
%            of zero enforces corresponding QOI constraint, i.e. allows no
%            relaxation.
%     wX    - A obj.Variables.Length-by-2 matrix of weights to the
%            relaxations on the
%            variables.
%     nInit - number of initial conditions for fmincon, sampled from SDP
%            result. In some cases, multiple initial conditions will
%            improve the quality of the returned solution.
%     opt   - a B2BDC option
% --------------------------------------------------------------------------
% Instead of specifying a weight matrix, setting wY and/or wX to the
% character strings 'perc','uwidth','unit', or 'null' automatically activates one
% of the standard weighting configurations (percentage change, uncertainty
% width, absolute change, or zero)



%Step 1: Process inputs
nVar = dsObj.Variables.Length; %number of variables
nExp = dsObj.DatasetUnits.Length; %number of experiments
yBnds = [[dsObj.DatasetUnits.Values.LowerBound]' [dsObj.DatasetUnits.Values.UpperBound]']; %[lb ub measurement];
xBnds = [[dsObj.Variables.Values.LowerBound]' [dsObj.Variables.Values.UpperBound]'];

if ischar(wY) && strcmpi(wY, 'perc') %percentage consistency measure
    wY = abs(yBnds);
elseif ischar(wY) && strcmpi(wY, 'uwidth') %symmetric uncertainy width
    wY = repmat((yBnds(:,2) - yBnds(:,1)),1,2);
elseif ischar(wY) && strcmpi(wY, 'unit')
    wY = ones(nExp,2);
elseif ischar(wY) && strcmpi(wY, 'null') || norm(wY) == 0;
    wY = zeros(nExp,2);
elseif size(wY,1) ~= nExp
    error('Improper QOI weight specification')
end
if ischar(wX) && strcmpi(wX, 'perc')
    wX = abs(xBnds);
elseif ischar(wX) && strcmpi(wX, 'uwidth') %symmetric uncertainty width
    wX = repmat((xBnds(:,2) - xBnds(:,1)),1,2);
elseif ischar(wX) && strcmpi(wX, 'unit')
    wX = ones(nVar,2);
elseif ischar(wX) && strcmpi(wX, 'null') || norm(wX)==0;
    wX = zeros(nVar,2);
elseif size(wX,1) ~= nVar
    error('Improper variable weight specification')
end

if nargin < 4
    nInit =1;
end

if nargin < 5
   opt = generateOpt;
end
units = dsObj.DatasetUnits.Values;
absErr = zeros(dsObj.Length,1);
if opt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats)
         absErr(j) = units(j).SurrogateModel.ErrorStats.absMax;
      end
   end
end

% Ignore noninfluential/redundant relaxations, i.e. those with zero weight.
% Used to build the problem with the minimum number of dimensions needed.

indYLb = find(wY(:,1)>0);
indYUb = find(wY(:,2)>0);
indXLb = find(wX(:,1)>0);
indXUb = find(wX(:,2)>0);

nYLb = numel(indYLb);
nYUb = numel(indYUb);
nXLb = numel(indXLb);
nXUb = numel(indXUb);

wyl = wY(indYLb,1);
wyu = wY(indYUb,2);
wxl = wX(indXLb,1);
wxu = wX(indXUb,2);

% Step 2: Build relevant objective and constraint matrices
dim = 1+nVar+nYLb+nYUb+nXLb+nXUb; %[1 x DL DU dl du], dimension of the problem
varList = {dsObj.Variables.Values.Name}';
[Cobj, dsuConstr, varConstr] = buildObjAndConstraintMatrices;
constrSet = [dsuConstr; varConstr];

% Step 3: Sedumi implementation of the optimization in primal form
% min c'x : Ax = b, x \in K
nConstr = numel(constrSet);
p = nConstr+dim^2;

b = [spalloc(nConstr,1,0); 1];
c = [spalloc(nConstr,1,0); Cobj(:)];

nzMax = 1+nConstr;
for kk = 1:nConstr
    temp = constrSet{kk}';
    nzMax = nzMax+nnz(temp);
end

A = spalloc(1+nConstr,p, nzMax);
A(1:nConstr, 1:nConstr) = speye(nConstr);
for kk = 1:nConstr
    temp = constrSet{kk}';
    A(kk,nConstr+1:p) = temp(:)';
end
A(1+nConstr,1+nConstr) = 1;

K.l = nConstr;
K.s = dim;
pars.fid = 0;
pars.eps = 0; 
[x,y, info] = sedumi(A,b,c,K,pars);
sdpObj = c'*x;
Z = mat(x(nConstr+1:end));
z = Z(2:end,1);
E = Z(2:end,2:end);

% Step 5: Finding an upper bound via fmincon
CovZ = E-z*z';
func = @(xx) objFn(xx,Cobj);
%ineqC = @(xx) constrFn(xx,dsuConstr);
ineqC = @(xx) constrFn(xx,constrSet);

zLB = -inf*ones(nVar+nYLb+nYUb+nXLb+nXUb,1);
zUB = inf*ones(nVar+nYLb+nYUb+nXLb+nXUb,1);

zSamp = repmat(z,1,nInit-1) + chol(CovZ+(1e-5)*eye(dim-1))*randn(numel(z),nInit-1);
%zSamp = lhsdesign(nInit-1, dim-1);
z0 = [z zSamp];

% options = optimset('Display', 'on','Algorithm','interior-point',...
%     'MaxFunEvals',10000,'GradObj','on') ;
% options = optimset('Display', 'on','Algorithm','interior-point',...
%     'MaxFunEvals',10000,'GradObj','on','DerivativeCheck','on', ...
%     'FinDiffType', 'central','GradConstr','on','Hessian','user-supplied', ...
%     'HessFcn',@(x,lambda)quadhess(x,lambda,Cobj,constrSet)) ;
% options = optimset('Display', 'off','Algorithm','interior-point',...
%     'TolFun',1e-10,'TolCon', 0,'TolX', 0,'GradObj','on','GradConstr','on',...
%     'Hessian','user-supplied','HessFcn', ...
%     @(x,lambda)quadhess(x,lambda,Cobj,dsuConstr)) ;
options = optimset('Display', 'on','Algorithm','interior-point',...
    'TolFun',1e-12,'TolCon', 1e-12,'TolX', 1e-12,'MaxIter', 1e5,...
    'MaxFunEvals',1e5,'GradObj','on','GradConstr','on',...
    'Hessian','user-supplied','HessFcn', ...
    @(x,lambda)quadhess(x,lambda,Cobj,constrSet)) ;
objFeas = inf;
for kk = 1:nInit
    [zTemp,objTemp,EXITFLAG, fopt] = fmincon(func,z0(:,kk),[],[],[],[],zLB,zUB, ineqC, options);
    
    if  objTemp < objFeas
        zFeas = zTemp;
        objFeas = objTemp;
    end    
end
xFeas = zFeas(1:nVar);

%Restore relaxations to original size
DeltaLb = zeros(nExp,1);
DeltaLb(indYLb) = zFeas(nVar+1:nVar+nYLb);
DeltaUb = zeros(nExp,1);
DeltaUb(indYUb) = zFeas(nVar+nYLb+1:nVar+nYLb+nYUb);
deltaLb = zeros(nVar,1);
deltaLb(indXLb) = zFeas(nVar+nYLb+nYUb+1:nVar+nYLb+nYUb+nXLb);
deltaUb = zeros(nVar,1);
deltaUb(indXUb) = zFeas(nVar+nYLb+nYUb+nXLb+1:nVar+nYLb+nYUb+nXLb+nXUb);
% The reported relaxations are relative to the given weights.
vcReport.Objective.SDPLowerBound = sdpObj;
vcReport.Objective.Feasible = objFeas;
vcReport.FeasiblePoint = xFeas;
vcReport.Relaxations.yLowerBound = DeltaLb;
vcReport.Relaxations.yUpperBound = DeltaUb;
vcReport.Relaxations.xLowerBound = deltaLb;
vcReport.Relaxations.xUpperBound = deltaUb;
vcReport.Weights.yLowerBound = wY(:,1);
vcReport.Weights.yUpperBound = wY(:,2);
vcReport.Weights.xLowerBound = wX(:,1);
vcReport.Weights.xUpperBound = wX(:,2);

%Check
if any(ineqC(zFeas)>0)
    warning('Fmincon failed to find a feasible relaxation.');
    violation = norm(ineqC(zFeas))
end
% Nested functions
    function [Cobj, dsuConstr,varConstr] = buildObjAndConstraintMatrices() 
        % This function can be sped up by rewritting with sparse() as
        % opposed to spalloc()
        Q_U = cell(nExp,1); %Upperbnd on QOI
        Q_L = cell(nExp,1); %Lowerbnd on QOI
        P_U = cell(nYUb,1); %Nonnegativity of upperbnd relaxation
        P_L = cell(nYLb,1); %Nonnegativity of lowerbnd relaxation
        
        for jj = 1:nExp
            dsUnit = dsObj.DatasetUnits.Values(jj);
            numActVar = dsUnit.VariableList.Length;
            M = expandModels(dsUnit);   %expands the jjth model from active variables to full variables.
            Aquad = M(2:nVar+1, 2:nVar+1);
            const = M(1,1);
            blin = M(2:nVar+1,1);
            nzMax = 1 + numActVar^2+ 2*numActVar + 2;
            U = yBnds(jj,2)+absErr(jj);
            L = yBnds(jj,1)-absErr(jj);
            
            Qltemp = spalloc(dim,dim, nzMax);
            Qltemp(1,2:nVar+1) = -blin';
            Qltemp(2:nVar+1,1) = -blin;
            Qltemp(1,1) = -const+L;
            Qltemp(2:nVar+1, 2:nVar+1) = -Aquad;
            
            Qutemp = spalloc(dim,dim, nzMax);
            Qutemp(1,2:nVar+1) = blin';
            Qutemp(2:nVar+1,1) = blin;
            Qutemp(1,1) = const-U;
            Qutemp(2:nVar+1, 2:nVar+1) = Aquad;
            
            ind = find(jj==indYLb);
            if ~isempty(ind)
                Qltemp(1,nVar+1+ind) = -0.5*wyl(ind);
                Qltemp(nVar+1+ind,1) = -0.5*wyl(ind);
            end
            
            ind = find(jj==indYUb);
            if ~isempty(ind)
                Qutemp(1,nVar+1+nYLb+ind) = -0.5*wyu(ind);
                Qutemp(nVar+1+nYLb+ind,1) = -0.5*wyu(ind);
            end
            
            Q_U{jj} = Qutemp;
            Q_L{jj} = Qltemp;
        end
        
        for kk =1:nYLb
            Pltemp = spalloc(dim,dim,2);
            Pltemp(1,nVar+1+kk) = -0.5;
            Pltemp(nVar+1+kk,1) = -0.5;
            P_L{kk} = Pltemp;
        end
        
        for kk = 1:nYUb
            Putemp = spalloc(dim,dim,2);
            Putemp(1,nVar+1+nYLb+kk) = -0.5;
            Putemp(nVar+1+nYLb+kk,1) = -0.5;
            P_U{kk} = Putemp;
        end
        
        Q_H = cell(nVar,1); %Bounds on variables (quadratic, includes relaxations: lb - relLB <= x <= ub +relUB)
        P_u = cell(nXUb,1); %Nonnegativity of upperbnd relaxation
        P_l = cell(nXLb,1); %Nonnegativity of lowerbnd relaxation

        for ii = 1:nVar  %bounds on X
            l = xBnds(ii,1);
            u = xBnds(ii,2);
            Jx = zeros(nVar,nVar);
  
            nzMaxH = 15;
            Htemp = spalloc(dim,dim,nzMaxH);
            Htemp(1,1+ii) = -(l+u)/2;
            Htemp(1+ii,1) = -(l+u)/2;
            Htemp(1,1) = l*u;
            Htemp(1+ii,1+ii) = 1;
            
            ind1 = find(ii==indXLb);
            if ~isempty(ind1)
                
                Htemp(1,1+nVar+nYLb+nYUb+ind1) = -0.5*wxl(ind1)*u;
                Htemp(1+nVar+nYLb+nYUb+ind1,1) = -0.5*wxl(ind1)*u;
                Htemp(1+nVar+nYLb+nYUb+ind1,1+ii) = 0.5*wxl(ind1);
                Htemp(1+ii,1+nVar+nYLb+nYUb+ind1) = 0.5*wxl(ind1);
                
                Phltemp = spalloc(dim,dim,2);
                Phltemp(1,1+nVar+nYLb+nYUb+ind1) = -0.5;
                Phltemp(1+nVar+nYLb+nYUb+ind1,1) = -0.5;
                P_l{ind1} = Phltemp;
            end
            ind2 = find(ii==indXUb);
            if ~isempty(ind2)
                
                Htemp(1,1+nVar+nYLb+nYUb+nXLb+ind2) = 0.5*wxu(ind2)*l;
                Htemp(1+nVar+nYLb+nYUb+nXLb+ind2,1) = 0.5*wxu(ind2)*l;
                Htemp(1+nVar+nYLb+nYUb+nXLb+ind2,1+ii) = -0.5*wxu(ind2);
                Htemp(1+ii,1+nVar+nYLb+nYUb+nXLb+ind2) = -0.5*wxu(ind2);
                
                Phutemp = spalloc(dim,dim,2);
                Phutemp(1,1+nVar+nYLb+nYUb+nXLb+ind2) = -0.5;
                Phutemp(1+nVar+nYLb+nYUb+nXLb+ind2,1) = -0.5;
                P_u{ind2} = Phutemp;
            end
            
            if ~(isempty(ind1) || isempty(ind2))
                
                Htemp(1+nVar+nYLb+nYUb+ind1, ...
                    1+nVar+nYLb+nYUb+nXLb+ind2) = wxl(ind1)*wxu(ind2);
                Htemp(1+nVar+nYLb+nYUb+nXLb+ind2,...
                    1+nVar+nYLb+nYUb+ind1) = wxl(ind1)*wxu(ind2);
            end
            Q_H{ii} = Htemp;
        end
        
        dsuConstr = [Q_U; Q_L; P_U; P_L];
        varConstr = [Q_H; P_u; P_l];
        nzMax = 2*(nYLb+nYUb+nXLb+nXUb);
        Cobj = spalloc(dim,dim, nzMax);
        Cobj(1,nVar+2:dim) = 0.5*ones(1,nYLb+nYUb+nXLb+nXUb);
        Cobj(nVar+2:dim, 1) = 0.5*ones(1,nYLb+nYUb+nXLb+nXUb);
    end


    function [Q] = expandModels(dsUnit)
        %Expands coefficient matrix from active variableList to full variableList
        %Input dsUnit is a datasetUnit
        activeVar = {dsUnit.VariableList.Values.Name}';
        nzMax = size(activeVar,1)^2;
        Q = spalloc(nVar+1,nVar+1,nzMax);
        [~,~,index] = intersect(activeVar,varList,'stable');
        index1 = [1; index+1];
        Q(index1, index1) = dsUnit.SurrogateModel.CoefMatrix; 
    end

end


%Local functions (Fmincon)

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

