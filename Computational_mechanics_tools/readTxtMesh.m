function [X, T] = readTxtMesh(filename)

fileID = fopen(filename,'r');

%% Read general data
% 
numComp=4;
anArray = fscanf(fileID, '%d', [1,numComp]);

numNodes=anArray(1);
numElems=anArray(2);
spaceDim=anArray(3);
numElemNodes=anArray(4);

fprintf('%6d%6d%6d%6d\n',...
        numNodes, numElems, spaceDim, numElemNodes);

%% Read coordinates
% 
X=zeros(numNodes,spaceDim);

numComp=spaceDim+1;
for i=1:numNodes
    anArray=fscanf(fileID, '%f', [1,numComp]);
    node=int64(anArray(1));
    X(node,:)=anArray(2:numComp);
end

%% Write connectivities
% 
T=zeros(numElems,numElemNodes);

numComp=numElemNodes+1;
for i=1:numElems
    anArray = fscanf(fileID, '%d', [1,numComp]);
    elem=anArray(1);
    T(elem,:)=anArray(2:numComp);
end

fclose(fileID);
