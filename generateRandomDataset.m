function ds = generateRandomDataset(nQOI,nVar,actVar,seedID)
% DS = GENERATERANDOMDATASET(NQOI,NVAR,SEEDID) generates a random dataset with
% nQOI random targets at variable dimension nVar, with random seed number
% seedID. The variables are assumed to have uncertainty range [-1,1].

%  Created: June 28, 2018     Wenyu Li

if nargin > 3
   rng(seedID);
end
if nargin < 3 || isempty(actVar)
   actVar = nVar;
end
tVar = generateVar([],repmat([-1,1],nVar,1));
ds = generateDataset('Random dataset');
for i = 1:nQOI
   tmpID = randperm(nVar,actVar);
   Q = rand(actVar+1);
   Q = Q+Q';
   tModel = generateModel(Q,tVar.makeSubset(tmpID));
   LB = Q(1,1)-rand(1);
   UB = Q(1,1)+rand(1);
   tUnit = generateDSunit(['Random QOI ' num2str(i)],tModel,[LB UB]);
   ds.addDSunit(tUnit);
end
