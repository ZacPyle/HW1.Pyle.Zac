% This is a static force analysis of a 4th order Michell Truss. 
% This script a matrix containing xyz coordinates of the fixed and 
% free nodes, the connectivity matrix describing the connection between 
% these nodes, and applied loads. For all matrices with coordinates (such
% as P, Q, C the top row is x coordinates, 2nd row is y coordinates, and
% each column is a node/member
% 
% This truss is designed so that all members in compression are 
% bars and all members in tension are strings. 
% 
% Created: 10/30/19
% Author : Zac Pyle
% PID    : A12601746

clear all;
close all;
clc;

% ---------------------------------------------------
%   Basic Variables
% ---------------------------------------------------

Phi   = pi/16;             % angle between rl's radiating from origin
Beta  = pi/6;              % Angle between rl and pl
q     = 10;                % Number of free nodes (nodes between members)
p     = 5;                 % Number of fixed nodes (nodes on inner circle)
b     = 10;                % Number of bars (should always be in compression
s     = 10;                % Number of strings (should always be in tension
Theta = -pi/2:.1:6*pi/12;   % vector of theta values that will be useful for creating my circles
Order = 4;
dim   = 2;
n     = q + p;
m     = b + s;
Cq    = zeros(m, q);
Cp    = zeros(m, p);


% ---------------------------------------------------
%   Define properties of the truss
% ---------------------------------------------------

a = sin(Beta)/sin(Beta + Phi);
c = sin(Phi)/sin(Beta + Phi);

% Anchor circle of truss
% r4x = r4*cos(Theta);
% r4y = r4*sin(Theta);

eval("r" + Order + " = 1;");   % The radius of the anchor circle (circle the fixed nodes of the truss will be mounted on)
eval("r" + Order + "x = r" + Order + " * cos(Theta);");
eval("r" + Order + "y = r" + Order + " * sin(Theta);");

radMat = [eval("r" + Order)];

% ---------------------------------------------------
%   Create the Truss Outline
% ---------------------------------------------------

for i = Order:-1:1
    % Create radii
    temp1 = "r" + (i-1) + " = " + "r" + i + " / a;";    
    eval(temp1);  
    
    % Create xy coordinates for radii
    temp2 = "r" + (i-1) + "x = r" + (i-1) + " * cos(Theta);";
    temp3 = "r" + (i-1) + "y = r" + (i-1) + " * sin(Theta);";
    eval(temp2);
    eval(temp3);
end
    
% Gather xy coordinates of circles together into matrices and radii
for i = 0:Order
   xCircles(i+1,:) = eval( "r" + i + "x");
   yCircles(i+1,:) = eval( "r" + i + "y");
   
   radMat = [eval("r" + i); radMat];
end

% ---------------------------------------------------
%   Create P, Q, C matrices
% ---------------------------------------------------

% Create P Matrix (for fixed nodes on anchor circle)
% This lists fixed nodes starting from the bottom, moving CCW
% P(:,1) = [0;
%           0];
P(:,1) = [r4 * cos(-4*Phi);
          r4 * sin(-4*Phi)];
P(:,2) = [r4 * cos(-2*Phi);
          r4 * sin(-2*Phi)];
P(:,3) = [r4;
          0];
P(:,4) = [r4 * cos(2*Phi);
          r4 * sin(2*Phi)];
P(:,5) = [r4 * cos(4*Phi);
          r4 * sin(4*Phi)];
      
% Create Q Matrix (for free nodes on anchor circle)
% This matrix lists the nodes by arm, not by radius. 
% So each node in an arm is listed before switching 
% to the next arm.
arm = 0;
for j = 1:Order    
    for i = 1:(Order - j + 1)
    Q(1,arm + i) = radMat(i+1) * cos( (-4+((j-1)*2) + i)*Phi);
    Q(2,arm + i) = radMat(i+1) * sin( (-4+((j-1)*2) + i)*Phi);
    end
    arm = arm + (Order - (j - 1));
end

% Create the C matrix. Do this by first creating Cp
% and then creating Cq.
count = 1;
for i = 1:Order
    Cp(count, i) = -1;
    count = count + (Order-i+1);
end
for i = 1:Order
   Cp(count, end-i+1) = -1;
   count = count + (Order-i+1);
end

% for i = 1:p-1
%    Cp(i,1)    = -1;
%    Cp(i, i+1) = 1;
% end
% 
% count = p;
% for i = 1:Order
%    Cp(count, i+1) = -1;
%    count = count + (Order - i+1);
% end


% ---------------------------------------------------
%   Plot Truss for visualization purposes
% ---------------------------------------------------

figure(1);
hold on;
grid on;
for i = 0:Order
    if i == Order
        plot(xCircles(i+1,:), yCircles(i+1,:), 'k', 'LineWidth', 2);
    end
   plot(xCircles(i+1,:), yCircles(i+1,:), 'k'); 
end

for i = -4:4
    plot( [0 r0*cos(Phi*i)], [0 r0*sin(Phi*i)], 'k');
end

plot(P(1,:), P(2,:), 'mo');
plot(Q(1,:), Q(2,:), 'g*');
