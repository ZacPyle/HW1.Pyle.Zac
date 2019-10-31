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
Theta = -pi/2:.1:6*pi/12;  % vector of theta values that will be useful for creating my circles
Order = 4;                 % Order of the michell truss
dim   = 2;                 % 2D or 3D analysis
n     = q + p;             % Total number of nodes
m     = b + s;             % Total number of members
Cq    = zeros(m, q);       % Initialize various matrices for use later
Cp    = zeros(m, p);
P     = zeros(dim, p);
Q     = zeros(dim, q);
U     = zeros(dim, q);


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

% Create Bars portion of Cq
count = 1;
for i = 1:Order
    Cp(count, i) = -1;
    count = count + (Order-i+1);
end
for i = 1:Order
   Cp(count, end-i+1) = -1;
   count = count + (Order-i+1);
end

% Creat Cq
for i = 1:q
   Cq(i,i) = 1;
   if all(i ~= [4 7 9])
       Cq(i+1, i) = -1;
   end
end

% Just finish the strings of Cq by hand, too lazy to figure out a modular
% way to do it
Cq(11, 10) = 1;
Cq(12, 9)  = 1;
Cq(12, 10) = -1;
Cq(13, 7)  = 1;
Cq(13, 9)  = -1;
Cq(14, 7)  = -1;
Cq(14, 4)  = 1;
Cq(15, 8)  = 1;
Cq(16, 6)  = 1;
Cq(16, 8)  = -1;
Cq(17, 3)  = 1;
Cq(17, 6)  = -1;
Cq(18, 5)  = 1;
Cq(19, 2)  = 1;
Cq(19, 5)  = -1;
Cq(20, 1)  = 1;

% now combine Cq and Cp
C = [Cq Cp];


% -------------------------------------------------------------------
%   Create Applied Loads
% -------------------------------------------------------------------

% U matrix has vector forces applied to each free node, so it must
% be a dim x q matrix. When the loads are all pointing straight down
% that simulates the truss being mounted as seen in the graph. When
% loads are all to the right that simulates the truss mounted pointing
% straight down (rotated -90 degrees from as seen in the graph

% Only apply a load straight down at the tip of the truss (free node 4)
U1      = U;
U1(2,4) = -2;

% Only apply a load to the right at the tip of the truss (free node 4)
U2      = U;
U2(1,4) = 2;

% apply a unit load to all free nodes straight down
U3      = U;
U3(2,:) = -1;

% apply a unit load to all free nodes to the right
U4      = U;
U4(1,:) = 1;

%
% Now call tensegrity_Statics
%
[c_bars,t_strings,V] = tensegrity_statics(b,s,q,p,dim,Q,P,C,U1);
tensegrity_plot(Q,P,C,b,s,U1,V,true,2.0); 

% ---------------------------------------------------
%   Plot Truss for visualization purposes
% ---------------------------------------------------

% figure(1);
% hold on;
% grid on;
% for i = 0:Order
%     if i == Order
%         plot(xCircles(i+1,:), yCircles(i+1,:), 'k', 'LineWidth', 2);
%     end
%    plot(xCircles(i+1,:), yCircles(i+1,:), 'k'); 
% end
% 
% for i = -4:4
%     plot( [0 r0*cos(Phi*i)], [0 r0*sin(Phi*i)], 'k');
% end
% 
% 
% 
% plot(P(1,:), P(2,:), 'mo');
% plot(Q(1,:), Q(2,:), 'g*');
