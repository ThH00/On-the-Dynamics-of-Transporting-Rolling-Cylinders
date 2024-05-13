function [I_xbar,II_xbar] = get_config(R, h, theta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% To find a configuration where two cylinders are touching
% 1. Set cylinder I in a fixed position and orientation.
% 2. Set cylinder II in a fixed orientation.
% 3. Position cylinder II such that it is away from cylinder I so that the
% two cylinders do not intersect.
% 4. Plot the two cylinders

% 1. Set cylinder I in a fixed position and orientation.
I_psi = 0;
I_theta = theta;
I_phi = 0;

I_x1 = 0;
I_x2 = 0;
I_x3 = 0;

I_e2pp = [-cos(I_theta)*sin(I_psi);
    cos(I_theta)*cos(I_psi);
    sin(I_theta)];

I_e3 = [sin(I_theta)*sin(I_psi);
    -sin(I_theta)*cos(I_psi);
    cos(I_theta)];

I_xbar = [I_x1;I_x2;I_x3]+h/2*I_e3+R*I_e2pp;

% 2. Set cylinder II in a fixed orientation.
% 3. Position cylinder II such that it is away from cylinder I so that the
% two cylinders do not intersect.
II_psi = pi;
II_theta = theta;
II_phi = 0;

II_x1 = 3*R;
II_x2 = 0.5*R;
II_x3 = 0;

II_e2pp = [-cos(II_theta)*sin(II_psi);
    cos(II_theta)*cos(II_psi);
    sin(II_theta)];

II_e3 = [sin(II_theta)*sin(II_psi);
    -sin(II_theta)*cos(II_psi);
    cos(II_theta)];

II_x = [II_x1;II_x2;II_x3];
II_xbar = II_x+h/2*II_e3+R*II_e2pp;

% 4. Plot the two cylinders
figure()
hold on
subplot(1,2,1)
view(45,45)
hold on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
plot3([-1 1],[0 0],[0 0])
plot3([0,0],[-1 1],[0 0])
plot3([0,0],[0,0],[0 2])
[disktI,diskbI,cylinder_axisI,sidesI,pointAI,pointPI] = ...
    plot_cylinder(R,h,I_xbar(1),I_xbar(2),I_xbar(3),...
    I_psi,I_theta,I_phi,'green');
[disktII,diskbII,cylinder_axisII,sidesII,pointAII,pointPII] = ...
    plot_cylinder(R,h,II_xbar(1),II_xbar(2),II_xbar(3),...
    II_psi,II_theta,II_phi,'blue');

% 5. Find min distance
A = [1 -dot(II_e3,I_e3); -dot(II_e3,I_e3) 1];
B = [dot(II_xbar-I_xbar,I_e3); dot(I_xbar-II_xbar,II_e3)];

x = A\B;
m = x(1);
n = x(2);

% 6. Find and plot contact point
% point A on axis
I_xA = I_xbar+m*I_e3;
plot3(I_xA(1),I_xA(2),I_xA(3),'o','color','red')
II_xA = II_xbar+n*II_e3;
plot3(II_xA(1),II_xA(2),II_xA(3),'*','color','red')
% point B on rim
v = I_xA-II_xA;
u = v/norm(v);
I_xB = I_xbar+m*I_e3-R*u;
plot3(I_xB(1),I_xB(2),I_xB(3),'o','color','red')
II_xB = II_xbar+n*II_e3+R*u;
plot3(II_xB(1),II_xB(2),II_xB(3),'*','color','red')
% connect the two lines
plot3([I_xA(1) II_xA(1)],[I_xA(2) II_xA(2)],[I_xA(3) II_xA(3)],...
    'linewidth',2,'color','black')

% 7. contact test
t1 = m<h/2 && m>-h/2;
t2 = n<h/2 && n>-h/2;

% 8. Translation second cylinder
min_dist = norm(v);
% translation distance
t = min_dist-2*R;
II_x = II_x+t*u;
II_xbar = II_x+h/2*II_e3+R*II_e2pp;

% 9. Recalculate min distance
A = [1 -dot(II_e3,I_e3); -dot(II_e3,I_e3) 1];
B = [dot(II_xbar-I_xbar,I_e3); dot(I_xbar-II_xbar,II_e3)];

x = A\B;
m = x(1);
n = x(2);

% 10. Replot cylinders
subplot(1,2,2)
hold on
view(45,45)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
plot3([-1 1],[0 0],[0 0],'color','blue')
plot3([0,0],[-1 1],[0 0],'color','green')
plot3([0,0],[0,0],[0 2])
[disktI,diskbI,cylinder_axisI,sidesI,pointAI,pointPI] = ...
    plot_cylinder(R,h,I_xbar(1),I_xbar(2),I_xbar(3),...
    I_psi,I_theta,I_phi,'green');
[disktII,diskbII,cylinder_axisII,sidesII,pointAII,pointPII] = ...
    plot_cylinder(R,h,II_xbar(1),II_xbar(2),II_xbar(3),...
    II_psi,II_theta,II_phi,'blue');

% 6. Find and plot contact point
% point A on axis
I_xA = I_xbar+m*I_e3;
plot3(I_xA(1),I_xA(2),I_xA(3),'o','color','red')
II_xA = II_xbar+n*II_e3;
plot3(II_xA(1),II_xA(2),II_xA(3),'*','color','red')
% point B on rim
v = I_xA-II_xA;
min_dist = norm(v);
u = v/norm(v);
I_xB = I_xbar+m*I_e3-R*u;
plot3(I_xB(1),I_xB(2),I_xB(3),'o','color','red')
II_xB = II_xbar+n*II_e3+R*u;
plot3(II_xB(1),II_xB(2),II_xB(3),'*','color','red')
% connect the two lines
plot3([I_xA(1) II_xA(1)],[I_xA(2) II_xA(2)],[I_xA(3) II_xA(3)],...
    'linewidth',2,'color','black')

end

