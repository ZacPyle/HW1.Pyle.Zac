function [c_bars,t_strings,V]=tensegrity_statics(b, s, q, p, dim, Q, P, C, U, tau_min, w, extra_constraints)
% Static equilibrium analysis of tensegrity structures with fixed nodes
% and external forces, assuming the tension in all tendons are measureable.
% If the conditions of static equilibrium are fully determined, solve for it.
% Otherwise (if the conditions of static equilibrium are underdetermined):
% a) determine whether or not the system is pretensionable
%    (with tau>=tau_min in all strings for u=zero). Then:
% b1) if system is pretensionable, find tensioning of the strings with
%    tau>=tau_min in all strings for applied U while minimizing w'*x_tilde, or
% b2) if the system is not pretensionable, find tensioning of the strings
%    that maximizes the minimum tension.
% Finally, compute the reaction forces V.
% Code by Thomas Bewley (UCSD), written in Summer 2019 as a JPL Faculty Fellow.
%
% For further explanation, see section 2.1 of B2019, which this code follows closely:
% Bewley (2019) Stabilization of low-altitude balloon systems, Part 2:
% riggings with multiple taut ground tethers, analyzed as tensegrity systems;
%
% INPUTS:
% b=number of bars,
% s=number of strings
% q=number of free nodes,
% p=number of fixed nodes
% dim=dimension of system (2 or 3)
% Q=Q_(dim x q)=locations of free nodes
% P=P_(dim x p)=locations of fixed nodes
% C=C_(m x n)  =connectivity matrix
% U=U_(dim x q)=applied external force vector at each free node
% OUTPUTS:
% c_bars     =magnitude of the compression force in each bar (neg for tension)
% t_strings  =magnitude of the tension force in each string (neg for compression)
% V_(dim x p)=reaction force vector at each fixed node
% extra_constraints=optional argument, used if pulleys equilibrate 2 or more tether tensions

m = b + s; % number of members
n = q + p; % number of nodes

% set some default values
if nargin < 10 
    tau_min = 0.1;
    w       = ones(m, 1);
    extra_constraints = "";
elseif nargin < 11
    w = ones(m, 1);
    extra_constraints = "";
elseif nargin < 12
    extra_constraints = "";
end 
      
% partition Connectivity matrix C
CQ = C(:,1:q);
CP = C(:,q+(1:p));
N  = [Q P];                         % assemble Node matrix N

% compute and partition Member matrix M, based on C and N
M = N*C'; 
B = M(:,1:b); 
S = M(:,b+(1:s));

for i = 1:m 
    ell(i) = norm(M(:,i));            % lengths of members
    D(:,i) = M(:,i)/ell(i);           % Normalized members(Direction) matrix
end

x   = sym('x',[1 m]);                 % set up (3a) in B2019 symbolically
X   = diag(x);        
sys = D*X*CQ-U;                       % note that sys=sys_(dim x q). We seek diagonal X s.t. sys=0.

% set up x1 to xm as symbolic variables, convert (3a) to (3b) [i.e., Ase*x=u]
for i = 1:m 
    exp = "syms x" + i; 
    eval(exp);                        % set up a symbolic equationsToMatrix command in SYS
end

SYS = 'equationsToMatrix(['; 
for i = 1:dim 
    for j = 1:q 
        SYS = SYS + "sys(" + i + "," + j + ") == 0";
       if i < dim || j < q
           SYS = SYS + ",";
       end
    end 
end

SYS = SYS + extra_constraints + "],[";

for i = 1:m
    SYS = SYS + "x" + i; 
    if i < m 
        SYS = SYS + ","; 
    end
end

SYS=SYS+"])";

[Ase, u] = eval(SYS);         % finally, execute the symbolic equationsToMatrix command assembled above
Ase = eval(Ase); 
u   = eval(u);                % convert Ase and u to a regular matrix and vector

mhat = dim*q; 
nhat = m;                     % note that Ase=Ase_(mhat x nhat)

[U_,S_,V_] = svd(Ase);
s_ = diag(S_);                % compute SVD of Ase (see sec 2.1.1 of B2019)

for i = 1:length(s_)
    if s_(i) > 1e-8
        r = i; 
    end
end

r 

DOF = nhat-r;                 % determine rank and DOF of solution

U_underbar = U_(:,1:r); 
U_overbar  = U_(:,r+1:mhat); 
S_underbar = S_(1:r,1:r);     % partition SVD
V_underbar = V_(:,1:r); 
V_overbar  = V_(:,r+1:nhat); 

% determine pseudoinverse of Ase
Ase_plus   = V_underbar * inv(S_underbar) * U_underbar';  

if r < mhat
    disp('Warning: Ase is potentially inconsistent, implying the presence of soft modes,')
    disp('or instability!  More strings or fixed points should fix the problem.')
else
    disp("Ase is not potentially inconsistent. Good.")
end

disp(' ')

disp('Bar compressions and string tensions with loads as specified,')
disp('least squares solution (i.e., NO pretensioning):')

x        = Ase_plus * u; 
residual = norm(Ase*x-u);  % Solve Ase*x=u using pseudoinverse, check residual (error in your solution)

if r < mhat
   if residual < 1e-7
       disp("u in column space of Ase, so at least one solution exists, with residual " + residual + ".");
   else
       disp("u not in column space of Ase.  No solutions without deformation (residual " + residual + ").]");
       return
   end
end

[c_bars, t_strings] = extract_tension_and_compression(x,b,s);

options = optimoptions('linprog','Display','none');
if DOF > 0
    disp("Ase is underdetermined with "+DOF+" DOF. Checking now to see if system is pretensionable,");
    disp("with tension >= "+tau_min+" in all tethers for zero applied load."); disp(' ');
    
    % Start with x=0 solution, corresponding to u=0. System is pretentionable iff a
    % solution to the LP (5a)-(5b) exists with all strings in tension.
    V_tilde = V_overbar(b+1:b+s,:);
    w       = ones(s,1); 
    w_tilde = -V_tilde' * w;
    w_bar   = [w_tilde; 
              -w_tilde];
    A_bar   = [-V_tilde V_tilde];      % Again, we are solving the LP in the u=0 case here.
    b_bar   = -tau_min*ones(s,1);      % (We will redo this below for the specified u.)
    
    [c_bar, FVAL, exitflag] = linprog(-w_bar,A_bar,b_bar,[],[],[],[],options);
    
    if exitflag==1
        c = c_bar(1:DOF) - c_bar(DOF+1:2*DOF);
        disp('Result with external load ZERO, pretensioned with given tau_min');
        disp('while minimizing the L1 norm of the tensions:');
        
        x0 = V_overbar * c; 
        [c_bars, t_strings] = extract_tension_and_compression(x0,b,s);
        
        pretensionable = true;
        disp('Pretensionable!');  % Section 2.1.4 of B2019
        disp('Results with external forces u as specified and tensioned with given tau_min')
        disp('while minimizing the L1 norm of the tensions.');
        w       = ones(s,1); 
        w_tilde = -V_tilde' * w; 
        w_bar   = [w_tilde; 
                  -w_tilde];
        A_bar   = [-V_tilde V_tilde];                          % Set up the LP in (5a)-(5b)
        b_bar   = x(b+1:b+s) - tau_min * ones(s,1);            % to tension the strings
        c_bar   = linprog(-w_bar,A_bar,b_bar,[],[],[],[],options);
        c       = c_bar(1:DOF)-c_bar(DOF+1:2*DOF);
        x       = x + V_overbar*c; 
        
        [c_bars, t_strings]=extract_tension_and_compression(x,b,s);
    else
        pretensionable = false; disp('Not pretensionable!'); % Section 2.1.5 of B2019
        disp('Results with external forces u as specified and tensioned to maximize tau_min:')
        
        w_bar = [zeros(2*DOF,1); 1];
        A_bar = [-V_tilde V_tilde ones(s,1)];                % Set up the modified LP in (5c)
        b_bar = x(b+1:b+s);                                  % to tension the strings
        
        c_bar = linprog(-w_bar,A_bar,b_bar,[],[],[],[],options);
        c     = c_bar(1:DOF) - c_bar(DOF+1:2*DOF);
        x     = x+V_overbar * c; 
        
        [c_bars, t_strings]=extract_tension_and_compression(x,b,s);
    end
else
    disp("Ase is not underdetermined (thus, it is not tensionable). The above solution is unique.")
end
V = D * diag(x) * CP;  % Compute reaction forces at fixed points.
end % function tensegrity_statics.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c_bars, t_strings] = extract_tension_and_compression(x,b,s)
c_bars(1:b) = -x(1:b);

if min(c_bars)>=0 
    disp("No bars under tension.  Good.");
else 
    disp("Note: some bars not under compression.  Maybe replace them with strings?"), 
end

t_strings(1:s) = x(b+(1:s)); 
tau_min        = min(t_strings);

if min(tau_min) >0
    disp("The " + s + " strings are all under tension with tau_min=" + tau_min + ". Good.")
else 
    disp("Some strings not under tension. Needs different tensioning or external loads.") 
end 
disp(' ')
end % function extract_tension_and_compression
