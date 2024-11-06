%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to compute several quality measures for 2D meshes /quads & tris)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

%% Read imput mesh in dat format
% 

folder='meshes';
fileName='tris_optional';
fileType='txt';
inputFile=strcat(folder,'/',fileName,'.',fileType)
[X,T]=readTxtMesh(inputFile)

%% Plot initial mesh
%

elemType=0;
figureNumber=1;
PlotMesh(T,X,elemType,figureNumber)

%% Compute quality measure based on inradiuas and circimradius for triangles
% 

numElementNodes=size(T,2);

if numElementNodes==3

[Q, minQ, maxQ, meanQ] = qualityRadiusRatio(X,T);
[q,minq,maxq,meanq] = qualityShapeTris(X,T);

elseif numElementNodes==4

[q,minq,maxq,meanq] = qualityShapeQuads(X,T);
else
    error('Invalid number of element nodes')
end
