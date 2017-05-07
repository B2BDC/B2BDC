%% Arun and Andy, December 18, 2015

nX = 9;
nT = 12;
maxNameLength = 7;

asciiValues = [48:57 65:90 95 97:122];
nAV = numel(asciiValues);

varNames = cell(nX,1);
for i=1:nX
    nl = ceil((maxNameLength-2)*rand)+2;
    rp = randperm(nAV);
    varNames{i} = char(asciiValues(rp(1:nl)));
end

xFeas = 2.0*(rand(nX,1)-0.5);

MatsVars = struct('Matrix',cell(nT,1),'Vars',cell(nT,1),'YBounds',cell(nT,1));
for i=1:nT
    nVars = ceil((nX-1)*rand)+1;
    rp = randperm(nX);
    MatsVars(i).Vars = varNames(rp(1:nVars));
    M = randn(1+nVars);
    MatsVars(i).Matrix = M + M';
    [~,~,indexIntoX] = intersect(MatsVars(i).Vars, varNames, 'stable');
    xLocal = xFeas(indexIntoX);
    y = [1;xLocal]'*MatsVars(i).Matrix*[1;xLocal];
    MatsVars(i).YBounds = [y-(abs(y)+0.1)*0.6*rand y+(abs(y)+0.1)*0.6*rand];
end
save Friday18December MatsVars

