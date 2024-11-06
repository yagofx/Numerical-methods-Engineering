function [q, minq, maxq, meanq] = qualityShapeTris(X,T)

% Number of elements
n = length(T(:,1));

nn = 2; % space dimension 2 for triangles

%Initialization of variables to compute measures
W = [1, 0.5; 0, sqrt(3)/2];
q = zeros(n,1);

for i =1:n
     % Get the vertices of the i-th triangle
     v1 = X(T(i, 1), :);
     v2 = X(T(i, 2), :);
     v3 = X(T(i, 3), :);

     % Compute A matrix
     A = [v2(1)-v1(1), v3(1)-v1(1); v2(2)-v1(2), v3(2)-v1(2)];

     % Compute S matrix
     S = A*inv(W);

     % Compute sigma
     sigma = det(S);
     
     % Compute norm
     Frob = norm(S,'fro');

     % Compute shape quality measure
     q(i) = (nn*sigma^(2/nn))/(Frob^2);
end
minq = min(q);
maxq = max(q);
meanq = mean(q);

% Plot the histogram of Q
figure;
histogram(q, 'BinEdges', linspace(min(q), max(q), 10), 'DisplayStyle', 'bar');
title('Quality Shape Tris');
xlabel('Quality');
ylabel('NumElemes')

end

