%--------------------------------------------------------------------------
%
%       Homework 1 - Finite elements
%
%--------------------------------------------------------------------------

clear all,clc,close all

% Problem 1
% Define the symbolic variable
syms x u(x)

% Define the differential equation
ode = diff(u, x, 2) == -sin(x);

% Solve the differential equation with boundary conditions
u_sol = dsolve(ode, u(0) == 0, u(1) == 1);

% Display the solution
disp('The solution for u(x) is:')
disp(u_sol)

% Plot the solution
uu = [0; sin(0.5)-0.5*sin(1)+1/2;1];
figure()
fplot(u_sol, [0 1])
hold on
plot([0 0.5 1],uu)
xlabel('x','Interpreter','latex')
ylabel('u(x)','Interpreter','latex')
legend('Analytical solution', 'FEM solution')
title('Problem 1 validation','Interpreter','latex')
hold off
grid on



% Problem 2

% Define symbolic variables
syms x C1 C2 C3 C4

% Define constants
k = 4;    % Coefficient in the equation
Q1 = 10;  % Heat source term for x in [0, 4]

% Solve for region 1: x in [0, 4] with Q = 10
phi1_double_prime = -Q1 / k;
phi1_prime = int(phi1_double_prime, x) + C1;  % First integration
phi1 = int(phi1_prime, x) + C2;               % Second integration

% Solve for region 2: x in [4, 8] with Q = 0
phi2_double_prime = 0;
phi2_prime = int(phi2_double_prime, x) + C3;  % First integration
phi2 = int(phi2_prime, x) + C4;               % Second integration

% Apply boundary conditions
% 1. phi(0) = 0
eq1 = subs(phi1, x, 0) == 0;

% 2. Continuity at x = 4: phi1(4) = phi2(4)
eq2 = subs(phi1, x, 4) == subs(phi2, x, 4);

% 3. Continuity of phi' at x = 4: phi1'(4) = phi2'(4)
eq3 = subs(phi1_prime, x, 4) == subs(phi2_prime, x, 4);

% 4. No heat flux at x = 8: -4*phi2'(8) = 0 -> phi2'(8) = 0
eq4 = subs(phi2_prime, x, 8) == 0;

% Solve for constants C1, C2, C3, and C4
solution = solve([eq1, eq2, eq3, eq4], [C1, C2, C3, C4]);

% Substitute the constants back into phi1 and phi2
phi1_sol = simplify(subs(phi1, [C1, C2], [solution.C1, solution.C2]));
phi2_sol = simplify(subs(phi2, [C3, C4], [solution.C3, solution.C4]));

% Display the solutions for both regions
disp('Solution for phi in the region 0 <= x <= 4:')
disp(phi1_sol)
disp('Solution for phi in the region 4 < x <= 8:')
disp(phi2_sol)

% Plot the solution
% figure()
% fplot(phi1_sol, [0, 4], 'LineWidth', 1.5)
% hold on
% fplot(phi2_sol, [4, 8], 'LineWidth', 1.5)
% xlabel('x')
% ylabel('\phi(x)')
% title('Analytical Solution for 4\phi''''(x) + Q(x) = 0')
% legend('\phi(x) for 0 \leq x \leq 4', '\phi(x) for 4 < x \leq 8')
% grid on
% hold off


% Material properties
k = 4;

%Problem boundary conditions
phi_1 = 0;
q_end = 0;
Q = 10; %Heat source

% Discretization
l = 8;
nelem = 4;
nnodes = nelem + 1;
x=0:2:l;

% Nodal values
phi = zeros(nnodes,1);
phi(1) = phi_1;
le = l/nelem;

% Elemental stiffnes matrix
k_e = k/le.*[1 -1;-1 1];

% Global stiffnes matrix
k = zeros(nnodes);

for i=1:nelem
    k(i:i+1,i:i+1)=k_e+k(i:i+1,i:i+1);
end

%Reduced source vector
f=zeros(nnodes,1);
fq = Q*le/2;

heat_nodes = 1:(l/2/le+1);

for j=heat_nodes
    if j==heat_nodes(1) || j==heat_nodes(end)
    f(j)=fq+f(j);
    else
    f(j)=2*fq+f(j);
    end
end

%Add Neumaan boundary condition, heat flux
f(end) = f(end) -q_end;
f_r = f(2:end);

% Get reduced stiffnes matrix
K_r = k(2:end,2:end);

% Compute nodal teperatures
phi_r = K_r\f_r;
phi(2:end) = phi_r;

% Plot nodal temperatures
figure()
hold on
fplot(phi1_sol, [0, 4])
fplot(phi2_sol, [4, 8])
plot(x,phi)
xlabel('x','Interpreter','latex')
ylabel('\phi(x)','Interpreter','tex')
legend('Analytical solution','', 'FEM solution')
grid on
title('Problem 2 validation','Interpreter','latex')
hold off






