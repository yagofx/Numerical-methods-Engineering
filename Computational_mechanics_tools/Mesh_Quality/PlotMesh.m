function PlotMesh(T,X,elem,figureNumber,str,nonum)
% 
% PlotMesh(T,X,elem,str,nonum)
%   X:           nodal coordinates
%   T:           connectivities 
%   elem:        element type (0: quadrilateral, 1: triangle)
%   figureNumber the number of the figure
%   str:         linestyle, color and marker used in the plot (optional)
%   nonum        = 1 to show nodes' number(optional)



% Line style and color
if nargin == 4
    str1 = 'bx';
    str2 = 'b-';
else
    if str(1) == ':' | str(1) == '-'  
        str1 = 'bx'; 
        str2 = ['b' str];
    else
        str1 = [str(1) 'x'];
        str2 = str;
    end
end
 
nen = size(T,2);

if elem == 0 
    if nen <= 4
        NodeOrder = [1:nen,1];
    elseif nen == 9 
         NodeOrder = [1  5  2  6  3  7  4  8  1];
    elseif nen == 16
        NodeOrder = [1  5  6  2  7  8  3  9  10  4  11  12  1]; 
    elseif nen == 25
        NodeOrder = [1  5  6  7  2  8  9  10  3  11  12  13  4  14  15  16  1]; 
    end
elseif elem == 1 
    if nen <= 3
        NodeOrder = [1:nen,1];
    elseif nen == 4
        NodeOrder = [1  2  3  1];          
    elseif nen==6 || nen == 7
        NodeOrder = [1  4  2  5  3  6  1];      
    elseif nen == 10
        NodeOrder = [1  4  5  2  6  7  3  8  9  1]; 
    elseif nen == 15
        NodeOrder = [1  4  5  6  2  7  8  9  3  10  11 12  1]; 
    end
elseif elem == 2 ||elem == 3
    NodeOrder = 1;
end

figure(figureNumber)
clf

% Nodes
plot(X(:,1),X(:,2),str1,'Linewidth',2)
hold on
% Elements
for j = 1:size(T,1)
    plot(X(T(j,NodeOrder),1),X(T(j,NodeOrder),2),str2,'LineWidth',2)
end

% Nodes' number
if nargin==5
   if length(nonum)>1
      for I=1:length(nonum)
         ien=nonum(I);
         text(X(ien,1),X(ien,2),int2str(ien))
      end
   else
      for I=1:size(X,1)
         text(X(I,1),X(I,2),int2str(I))
      end
   end
end
axis('equal')    
axis('off') 
hold off 
