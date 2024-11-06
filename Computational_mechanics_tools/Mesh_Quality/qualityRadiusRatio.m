function [Q, minQ, maxQ, meanQ] = qualityRadiusRatio(X,T)

% Number of elements
n = length(T(:,1));

%Initialization of variables to compute measures
A = zeros(n,1);
R = zeros(n,1);
r = zeros(n,1);
Q = zeros(n,1);

alpha = 2; % Normalization paramters

for i = 1:n
    % Get the vertices of the i-th triangle
     v1 = X(T(i, 1), :);
     v2 = X(T(i, 2), :);
     v3 = X(T(i, 3), :);
        
     % Compute the edge lengths
     L1 = norm(v2 - v3);
     L2 = norm(v3 - v1);
     L3 = norm(v1 - v2);
        
     % Semi-perimeter
     Sp = (L1 + L2 + L3) / 2;
        
     % Area of the triangle
     A(i) = 0.5*abs(det([v2(1)-v1(1), v3(1)-v1(1);v2(2)-v1(2), v3(2)-v1(2)]));
        
     % Circumradius R
     R(i) = (L1*L2*L3)/(4*A(i));
        
     % Inradius r
     r(i) = A(i)/Sp;
        
     % Radius ratio quality measure
     Q(i) = alpha*r(i)/R(i);
end
minQ = min(Q);
maxQ = max(Q);
meanQ = mean(Q);

% Plot the histogram of Q
figure;
histogram(Q, 'BinEdges', linspace(min(Q), max(Q), 10), 'DisplayStyle', 'bar');
title('Quality Radius Ratio Q');
xlabel('Quality');
ylabel('NumElemes')

end

