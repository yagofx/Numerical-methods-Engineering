%--------------------------------------------------------------------------
%
%               Homework 2: Integration
%
%--------------------------------------------------------------------------

clear all, clc, close all

format long


% Function definition
g1 = @g1;
g2 = @g2;
g3 = @g3;
g4 = @g4;

% Integral bounds
a = -4;
b = 4;

% Number of quadrature points
n = 5; % Example: 5 points for higher accuracy
disp('-----------------------------------------')
disp(['Integration points: ', num2str(n)]);
disp('-----------------------------------------')

% Get Legendre-Gauss nodes (xi) and weights (wi)
[xi, wi] = legendre_gauss_nodes_weights(n);

% Map nodes to the interval [a, b]
x_mapped = (b-a)/2 * xi + (a+b)/2;
weights_mapped = (b-a)/2 * wi;

% Compute integrals using the Legendre-Gauss quadrature formula
I1 =0;
I2 =0;
I3 =0;
I4 =0;
for i =1:n
I1 = weights_mapped(i)*g1(x_mapped(i)) + I1;
I2 = weights_mapped(i)*g2(x_mapped(i)) + I2;
I3 = weights_mapped(i)*g3(x_mapped(i)) + I3;
I4 = weights_mapped(i)*g4(x_mapped(i)) + I4;
end

% Display results
disp(['Integral of g1: ', num2str(I1)]);
disp(['Integral of g2: ', num2str(I2)]);
disp(['Integral of g3: ', num2str(I3)]);
disp(['Integral of g4: ', num2str(I4)]);
disp('-----------------------------------------')
% Errors
I1_r = 0.043662222222213;
I2_r = 1.494267689296229;
I3_r = 12.162401581657875;
I4_r = 31.817025833333332;

Err_abs = zeros(4,1);
Err_rel = zeros(4,1);

Err_abs(1) = abs(I1-I1_r);
Err_abs(2) = abs(I2-I2_r);
Err_abs(3) = abs(I3-I3_r);
Err_abs(4) = abs(I4-I4_r);

Err_rel(1) = Err_abs(1)/I1_r;
Err_rel(2) = Err_abs(2)/I2_r;
Err_rel(3) = Err_abs(3)/I3_r;
Err_rel(4) = Err_abs(4)/I4_r;

disp(['Eror_abs of g1: ', num2str(Err_abs(1))]);
disp(['Eror_abs of g2: ', num2str(Err_abs(2))]);
disp(['Eror_abs of g3: ', num2str(Err_abs(3))]);
disp(['Eror_abs of g4: ', num2str(Err_abs(4))]);
disp('-----------------------------------------')
disp(['Eror_rel of g1: ', num2str(Err_rel(1))]);
disp(['Eror_rel of g2: ', num2str(Err_rel(2))]);
disp(['Eror_rel of g3: ', num2str(Err_rel(3))]);
disp(['Eror_rel of g4: ', num2str(Err_rel(4))]);


%--------------------------------------------------------------------------
% Function to compute Legendre-Gauss nodes and weights
%--------------------------------------------------------------------------
function [x, w] = legendre_gauss_nodes_weights(n)
    % Generate nodes and weights for Legendre-Gauss quadrature
    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2)); % Recurrence relation coefficients
    T = diag(beta, 1) + diag(beta, -1);        % Tridiagonal matrix
    [V, D] = eig(T);                          % Eigenvalue decomposition
    x = diag(D);                              % Nodes (roots of the polynomial)
    w = 2 * (V(1,:).^2)';                     % Weights
end


