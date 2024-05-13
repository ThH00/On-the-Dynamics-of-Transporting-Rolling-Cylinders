% symbolic constraint calculation of two touching cylinders with a
% constraint independence check and an application of the Frobenius
% integrability criteria.

%% Dimensions
syms h r real

%% Center of mass position and velocity coordinates
syms q [1 3] real
syms qdot [1 3] real

%% 3-1-3 Euler Angles and resulting basis vectors for both cylinders
% psi, theta, phi
syms nu [1 3] real

% rotation tensors
R1 = [cos(nu(1)) sin(nu(1)) 0; -sin(nu(1)) cos(nu(1)) 0; 0 0 1];
R2 = [1 0 0; 0 cos(nu(2)) sin(nu(2)); 0 -sin(nu(2)) cos(nu(2))];
R3 = [cos(nu(3)) sin(nu(3)) 0; -sin(nu(3)) cos(nu(3)) 0; 0 0 1];

% fixed basis vectors
E1 = [1;0;0];
E2 = [0;1;0];
E3 = [0;0;1];

% moving basis vectors
e1_p = R1'*E1;
e2_pp = (R2*R1)'*E2;
e3 = (R3*R2*R1)'*E3;

%% Angular Velocity
syms nudot [1 3] real
omega = nudot(1)*E3+nudot(2)*e1_p+nudot(3)*e3;

%% Frobenius Integrability Criteria parameters
R = 3;     % number of constraints
ndof = 6;  % number of degrees of freedom
W = sym(zeros(R,ndof));
Pi = sym(zeros(R,1));   % array of constraints

%% Constraint: Cylinder I rolls without slipping
% the position vector from the center of mass of cylinder I to its point of
% contact with the ground
pi1 = -h/2*e3-r*e2_pp;
v1 = qdot'+cross(omega,pi1);

Pi(1) = v1(1);
Pi(2) = v1(2);
Pi(3) = v1(3);
% Pi(3) = I_qdot3+pi_1(3);

%% Frobenius Integrability Criteria
U = [qdot nudot];

for i = 1:R
    for j = 1:ndof
        W(i,j) = diff(Pi(i),U(j));
    end
end

% rank_W = rank(W)

%% Constructing the S matrices
Q = [q nu];
S = sym(zeros(ndof,ndof,R));
for B = 1:R
    for L = 1:ndof
        for K = 1:ndof
            S(L,K,B) = diff(W(B,L),Q(K))-diff(W(B,K),Q(L));
        end
    end
end
% S = simplify(S);

%% Find null solutions of Wx=0 (x=U)
null_space_W = null(W);

% choosing distinct common pairs:
% combination of 2 out of size(null_space_W,2)
num_combinations = nchoosek(size(null_space_W,2),2);
combination_pairs = nchoosek(1:size(null_space_W,2),2); 

% For all disntince column pair (a,b) of null_space_W, we calculate 
I = sym(zeros(R,num_combinations));
for B = 1:R
    for j = 1:num_combinations
        a = null_space_W(:,combination_pairs(j,1));
        b = null_space_W(:,combination_pairs(j,2));
        I(B,j) = simplify(a'*S(:,:,B)*b);
    end
end

%% Numerical Experiment
r_val = 0.1;
h_val = 1.27;
theta_val = pi/10;
q3_val = r_val*sin(theta_val)+h_val/2*cos(theta_val);
vals_val = [0,-0.101120139798576,q3_val,0,theta_val,h_val,r_val];
vars_val = [q1, q2, q3, nu1, nu2, h, r];
W_val = double(subs(W,vars_val,vals_val));
rank_val = rank(W_val);
S_val = double(subs(S,vars_val,vals_val));


%% Find null solutions of Wx=0 (x=U)
null_space_W = null(W_val); % 13x4 sym

% choosing distinct common pairs:
% combination of 2 out of size(null_space_W,2)
num_combinations = nchoosek(size(null_space_W,2),2);
combination_pairs = nchoosek(1:size(null_space_W,2),2); 

% For all disntince column pair (a,b) of null_space_W, we calculate 
I = sym(zeros(R,num_combinations));
for B = 1:R
    for j = 1:num_combinations
        a = null_space_W(:,combination_pairs(j,1));
        b = null_space_W(:,combination_pairs(j,2));
        I(B,j) = a'*S_val(:,:,B)*b;
    end
end
I = double(I);