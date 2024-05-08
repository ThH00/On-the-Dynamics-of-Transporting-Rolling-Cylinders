function [diskt,diskb,cylinder_axis,sides,pointA,pointP,basis,xP] ...
    = plot_cylinder(x1,x2,x3,psi,theta,phi,color)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% global R h
R = 0.4;
h = 5.08;


%% basis vectors
% Rotation matrices
R1 = [cos(psi), sin(psi), 0;
      -sin(psi), cos(psi), 0;
      0, 0, 1];
  
R2 = [1, 0, 0;
      0, cos(theta), sin(theta);
      0, -sin(theta), cos(theta)];

R3 = [cos(phi), sin(phi), 0;
      -sin(phi), cos(phi), 0;
      0, 0, 1];

% Fixed basis
E1 = [1;0;0];
E2 = [0;1;0];
E3 = [0;0;1];

% {E1, E2, E3} components of {e1',e2',e3'}
e1p = R1'*E1;
e2p = R1'*E2;
e3p = R1'*E3;

% {E1, E2, E3} components of {e1'',e2'',e3''}
e1pp = (R2*R1)'*E1;
e2pp = (R2*R1)'*E2;
e3pp = (R2*R1)'*E3;

% {E1, E2, E3} components of {e1,e2,e3}
e1 = (R3*R2*R1)'*E1;
e2 = (R3*R2*R1)'*E2;
e3 = (R3*R2*R1)'*E3;

%% plotting a circle in the horizontal plane
angle = linspace(0,2*pi,40);
circ1 = R*cos(angle);
circ2 = R*sin(angle);

% top disk
xcirct = x1+circ1*e1(1)+circ2*e2(1)+h/2*e3(1);     %  m
ycirct = x2+circ1*e1(2)+circ2*e2(2)+h/2*e3(2);     %  m
zcirct = x3+circ1*e1(3)+circ2*e2(3)+h/2*e3(3);     %  m
% bottom disk
xcircb = x1+circ1*e1(1)+circ2*e2(1)-h/2*e3(1);     %  m
ycircb = x2+circ1*e1(2)+circ2*e2(2)-h/2*e3(2);     %  m
zcircb = x3+circ1*e1(3)+circ2*e2(3)-h/2*e3(3);     %  m
% cylinder sides
xcirc = [xcirct' xcircb'];
ycirc = [ycirct' ycircb'];
zcirc = [zcirct' zcircb'];

% top disk
diskt = patch('xdata',xcirct,'ydata',ycirct,'zdata',zcirct,...
    'facecolor',[0, 0, 1],'linewidth',2);
% bottom disk
diskb = patch('xdata',xcircb,'ydata',ycircb,'zdata',zcircb,...
    'facecolor',[0, 0, 0],'linewidth',2);
% cylinder axis
cylinder_axis = line('xdata',[x1+h/2*e3(1),x1-h/2*e3(1)],...
    'ydata',[x2+h/2*e3(2),x2-h/2*e3(2)],...
    'zdata',[x3+h/2*e3(3),x3-h/2*e3(3)],...
    'linewidth',2,'color','red');
% cylinder sides
for i = 1:length(angle)
    sides(i) = line('xdata',xcirc(i,:),'ydata',ycirc(i,:),...
        'zdata',zcirc(i,:),'linewidth',1,'color',color);
end
% fixed point on top rim of cylinder (A)
pointA = plot3(xcirct(1),ycirct(1),zcirct(1),'linewidth',2,'color','red','marker','*');
% contact point path (P)
xP = [x1;x2;x3]-h/2*e3+R*e2pp;
pointP = plot3(xP(1),xP(2),xP(3),'linewidth',5,'color','black','marker','.');
% {e1, e2, e3} basis vectors
basis(1) = quiver3(x1+h/2*e3(1), x2+h/2*e3(2), x3+h/2*e3(3), e1(1), e1(2), e1(3),'r','linewidth',2);
basis(2) = quiver3(x1+h/2*e3(1), x2+h/2*e3(2), x3+h/2*e3(3), e2(1), e2(2), e2(3),'g','linewidth',2);
basis(3) = quiver3(x1+h/2*e3(1), x2+h/2*e3(2), x3+h/2*e3(3), e3(1), e3(2), e3(3),'k','linewidth',2);


end

