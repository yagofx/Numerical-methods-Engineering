function [q, minq, maxq, meanq] = qualityShapeQuads(X,T)

% Number of elements
n = length(T(:,1));

nn = 2; % space dimension 2 for triangles
p = 2; % penalty for low-quality elements


%Initialization of variables to compute measures
sub_tri = [1, 3, 4;
    2, 3, 4;
    1 ,2 ,4;
    1, 2, 3];

WW = [1, 0; 0, 1];
q = zeros(n,1);

for i =1:n
     % Get the vertices of the i-th quadrilateral
     v = X(T(i, :), :);

     etap = 0;

     for j=1:4
        v1 = v(sub_tri(j, 1), :);
        v2 = v(sub_tri(j, 2), :);
        v3 = v(sub_tri(j, 3), :);
        
        % Compute A matrix
        A = [v2(1)-v1(1), v3(1)-v1(1); v2(2)-v1(2), v3(2)-v1(2)];

        % Compute S matrix
        S = A*WW;
        
        % Compute sigma
        sigma = det(S);
         
        % Compute norm
        Frob = norm(S,'fro');

        % Compute eta^p
        etap = ((Frob^2)/(nn*sigma^(2/nn)))^p+etap;
     end

     eta = (1/4*etap)^(1/p);

     q(i) = 1/eta;
end
minq = min(q);
maxq = max(q);
meanq = mean(q);

% Plot the histogram of Q
figure;
histogram(q, 'BinEdges', linspace(min(q), max(q), 10), 'DisplayStyle', 'bar');
title('Quality Shape Quads');
xlabel('Quality');
ylabel('NumElemes')

end

