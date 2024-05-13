% symbolic constraint calculation of two touching cylinders with a
% constraint independence check and an application of the Frobenius
% integrability criteria.

%% Dimensions
syms h r real

%% Center of mass position and velocity coordinates
syms I_q [1 3] real
syms II_q [1 3] real
syms I_qdot [1 3] real
syms II_qdot [1 3] real

%% 3-1-3 Euler Angles and resulting basis vectors for both cylinders
syms I_nu [1 3] real    % I_psi, I_theta, I_phi
syms II_nu [1 3] real   % II_psi, II_theta, II_phi

% rotation tensors
I_R1 = [cos(I_nu(1)) sin(I_nu(1)) 0; -sin(I_nu(1)) cos(I_nu(1)) 0; 0 0 1];
I_R2 = [1 0 0; 0 cos(I_nu(2)) sin(I_nu(2)); 0 -sin(I_nu(2)) cos(I_nu(2))];
I_R3 = [cos(I_nu(3)) sin(I_nu(3)) 0; -sin(I_nu(3)) cos(I_nu(3)) 0; 0 0 1];

II_R1 = [cos(II_nu(1)) sin(II_nu(1)) 0; -sin(II_nu(1)) cos(II_nu(1)) 0;...
    0 0 1];
II_R2 = [1 0 0; 0 cos(II_nu(2)) sin(II_nu(2));...
    0 -sin(II_nu(2)) cos(II_nu(2))];
II_R3 = [cos(II_nu(3)) sin(II_nu(3)) 0; -sin(II_nu(3)) cos(II_nu(3)) 0;...
    0 0 1];

% fixed basis vectors
E1 = [1;0;0];
E2 = [0;1;0];
E3 = [0;0;1];

% moving basis vectors
I_e1_p = I_R1'*E1;
I_e2_pp = (I_R2*I_R1)'*E2;
I_e3 = (I_R3*I_R2*I_R1)'*E3;
II_e1_p = II_R1'*E1;
II_e2_pp = (II_R2*II_R1)'*E2;
II_e3 = (II_R3*II_R2*II_R1)'*E3;

%% Angular Velocity
syms I_nudot [1 3] real
syms II_nudot [1 3] real
I_omega = I_nudot(1)*E3+I_nudot(2)*I_e1_p+I_nudot(3)*I_e3;
II_omega = II_nudot(1)*E3+II_nudot(2)*II_e1_p+II_nudot(3)*II_e3;

%% Frobenius Integrability Criteria parameters
R = 9;      % number of constraints
ndof = 12;  % number of degrees of freedom
W = sym(zeros(R,ndof));
Pi = sym(zeros(R,1));   % array of constraints

%% Constraint: Cylinder I rolls without slipping
% the position vector from the center of mass of cylinder I to its point of
% contact with the ground
pi_1 = -h/2*I_e3-r*I_e2_pp;
v1 = I_qdot'+cross(I_omega,pi_1);

Pi(1) = v1(1); %dot(v1,E1);
Pi(2) = v1(2); %dot(v1,E2);
Pi(3) = v1(3); %dot(v1,E3);

%% Constraint: Cylinder II rolls without slipping
% the position vector from the center of mass of cylinder II to its point
% of contact with the ground
pi_4 = -h/2*II_e3-r*II_e2_pp;
v4 = II_qdot'+cross(II_omega,pi_4);

Pi(7) = v4(1); %dot(v4,E1);
Pi(8) = v4(2); %dot(v4,E2);
Pi(9) = v4(3); %dot(v4,E3);

%% Constraint: Cylinders I and II roll without slip wrt each other
A = [1, -dot(I_e3,II_e3);
    -dot(I_e3,II_e3), 1];
B = [dot(II_q-I_q,I_e3);
    dot(I_q-II_q,II_e3)];
params = A\B;
I_d = params(1);
II_d = params(2);

I_x = I_q'+I_d*I_e3;
II_x = II_q'+II_d*II_e3;

v = II_x-I_x;
norm_v = norm(v);

I_u = v/norm_v;
II_u = -I_u;

% position from the center of mass of cylinder I to the point of contact
% between cylinders I and II
pi_2 = I_d*I_e3+r*I_u;
% position from the center of mass of cylinder II to the point of contact
% between cylinders I and II
pi_3 = II_d*II_e3+r*II_u;

v2 = I_qdot'+cross(I_omega,pi_2);
v3 = II_qdot'+cross(II_omega,pi_3);

v_slip = v2-v3;

Pi(4) = v_slip(1);
Pi(5) = v_slip(2);
Pi(6) = v_slip(3);

%% Frobenius Integrability Criteria
u = [I_qdot I_nudot II_qdot II_nudot];
for i = 1:R
    for j = 1:ndof
        W(i,j) = diff(Pi(i),u(j));
    end
end
% symvar(W)

%% Constructing the S matrices
q = [I_q I_nu II_q II_nu];
S = sym(zeros(ndof,ndof,R));
for B = 1:R
    for L = 1:ndof
        for K = 1:ndof
            S(L,K,B) = diff(W(B,L),q(K))-diff(W(B,K),q(L));
        end
    end
end

%% Numerical Experiment
r_val = 0.1;
h_val = 1.27;
I_theta_val = pi/10;
II_theta_val = pi/10;
I_q3_val = -r_val*sin(I_theta_val)-h_val/2*cos(I_theta_val);
II_q3_val = -r_val*sin(II_theta_val)-h_val/2*cos(II_theta_val);
vals = [0,-0.101120139798576,I_q3_val,0,I_theta_val,...
    2*r_val,-0.198879860201424,II_q3_val,pi,II_theta_val,h_val,r_val];
vars = [I_q1, I_q2, I_q3, I_nu1, I_nu2,...
    II_q1, II_q2, II_q3, II_nu1, II_nu2, h, r];
W_val = double(subs(W,vars,vals));
S_val = double(subs(S,vars,vals));
rank_W = rank(W); % answer is 9

%% Find null solutions of Wx=0 (x=U)
null_space_W = null(W_val);

% choosing distinct pairs:
% combination of 2 out of size(null_space_W,2)
num_combinations = nchoosek(size(null_space_W,2),2);
combination_pairs = nchoosek(1:size(null_space_W,2),2); 

% For all distinct column pair (a,b) of null_space_W, we calculate 
I_val = sym(zeros(R,num_combinations));
for B = 1:R
    for j = 1:num_combinations
        a_val = null_space_W(:,combination_pairs(j,1));
        b_val = null_space_W(:,combination_pairs(j,2));
        I_val(B,j) = double(a_val'*S_val(:,:,B)*b_val);
    end
end
I_val = double(I_val);