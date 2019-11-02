% This is a static force analysis of a 4-Bar Nonminimal Tensegrity Prism. 
% This script a matrix containing xyz coordinates of the fixed and 
% free nodes, the connectivity matrix describing the connection between 
% these nodes, and applied loads. For all matrices with coordinates (such
% as P, Q, C the top row is x coordinates, 2nd row is y coordinates, and
% each column is a node/member. I have chosed to write the node
% coordninates using scale factors x, y, and z in case I want to scale one
% or multiple dimensions of this prism.
% 
% This truss is designed so that all members in compression are 
% bars and all members in tension are strings. 
% 
% The nodes are labeled as follows:
% 
%             5------8
%            /|     /|
%           / |    / |
%          6------7  |
%          |  |   |  |
%          |  |   |  |
%          |  1---|--4
%          | /    | /
%          |/     |/
%          2------3
%
% With node 1 at the origin
% 
% Created: 10/30/19
% Author : Zac Pyle
% PID    : A12601746

clear all;
close all;
clc;

% --------------------------------------------------
%   Basic Variables
% --------------------------------------------------

q   = 8;              % Number of free nodes
p   = 0;              % Number of fixed nodes
b   = 4;              % Number of bars
s   = 24;             % Number of strings
x   = 1;              % scale factor for x values
y   = 1;              % Scale factor for y values
z   = 1;              % scale factor for z values
n   = q + p;          % Total number of nodes
m   = b + s;          % Total number of members
dim = 3;              % Dimensions this problem has (2D or 3D)
P   = zeros(1,p);     % Initializes matrices that will be used later
Q   = zeros(1,q); 
U   = zeros(dim, q);
Cq  = zeros(m, q);
Cp  = zeros(m, p);

% --------------------------------------------------
%   Create P and Q Matrices
% --------------------------------------------------

% Since p = 0 we don't need to do anything, just keep it as an empty
% matrix. It was initialized above in case in the future this does have
% fixes nodes

% Just write out the Q matrix, not going to make it modular
Q(1,2) = x;
Q(1,3) = x;
Q(2,3) = y;
Q(2,4) = y;
Q(3,5) = z;
Q(1,6) = x;
Q(3,6) = z;
Q(:,7) = [x; y; z];
Q(2,8) = y;
Q(3,8) = z;
%Q(:,9) = 0.5*[x; y; z];

% --------------------------------------------------
%   Create C matric
% --------------------------------------------------

% Again, since p = 0 I don't have to do anything to Cp. It's just
% initialized incase in the future somneone wants to add fixed nodes
% NOTE: This will be a bit janky, there's more adaptable ways to code this
% but I went for the fastest

% Create bars section for Cq
Cq(1,1) = -1;
Cq(1,7) = 1;
Cq(2,4) = -1;
Cq(2,6) = 1;
Cq(3,3) = -1;
Cq(3,5) = 1;
Cq(4,2) = -1;
Cq(4,8) = 1;

% Create strings section of Cq on the bottom face of prism
for i = 1:3
   Cq(b+i, i)   = -1;
   Cq(b+i, i+1) = 1;
end
Cq(8,1) = 1;
Cq(8,4) = -1;

% Create strings section of Cq on sides of prism
% for i = 1:4
%    Cq(b+4+i, i)   = -1;
%    Cq(b+4+i, i+4) = 1;
% end

% Create strings section of Cq on top of prism
for i = 1:3
    Cq(b+8+i, i+4) = -1;
    Cq(b+8+i, i+5) = 1;
end
Cq(16,5) = -1;
Cq(16,8) = 1;

% Create the cross section strings on the sides now
Cq(17,1) = -1;
Cq(17,6) = 1;
Cq(18,2) = -1;
Cq(18,5) = 1;
Cq(19,2) = -1;
Cq(19,7) = 1;
Cq(20,3) = -1;
Cq(20,6) = 1;
Cq(21,3) = -1;
Cq(21,8) = 1;
Cq(22,4) = -1;
Cq(22,7) = 1;
Cq(23,4) = -1;
Cq(23,5) = 1;
Cq(24,1) = -1;
Cq(24,8) = 1;
 
% Create cross section strings on top and bottom, starting with bottom
Cq(25,1) = -1;
Cq(25,3) = 1;
Cq(26,2) = -1;
Cq(26,4) = 1;
Cq(27,5) = -1;
Cq(27,7) = 1;
Cq(28,6) = -1;
Cq(28,8) = 1;

C = [Cq Cp];

% -------------------------------------------------------------
%   Create U matrix (Applied loading matrix)
% -------------------------------------------------------------

% Loading for prism under compression
U1 = U;
U1(3,1:4) = 1;
U1(3,5:8) = -1;

% Loading for prism under tension
U2 = U;
U2(3,1:4) = -1;
U2(3,5:8) = 1;

% Loading on two opposite corners of prism. Vectors of unit magnitude,
% pointing towards center of cube
U3 = U;
U3(:,1) = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
U3(:,2) = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
U3(:,8) = [-1/sqrt(3); -1/sqrt(3); -1/sqrt(3)];
U3(:,7) = [-1/sqrt(3); -1/sqrt(3); -1/sqrt(3)];

% -------------------------------------------------------------
%   Now call tensegrity_Statics
% -------------------------------------------------------------
[c_bars,t_strings,V] = tensegrity_statics(b,s,q,p,dim,Q,P,C,U2);
tensegrity_plot(Q,P,C,b,s,U2,V,true,2.0); 

% --------------------------------------------------
%   Plot for Visualization
% --------------------------------------------------

% figure(1);
% hold on;
% grid on;
% 
% plot3(Q(1,:), Q(2,:), Q(3,:), 'go', 'MarkerSize', 12);


