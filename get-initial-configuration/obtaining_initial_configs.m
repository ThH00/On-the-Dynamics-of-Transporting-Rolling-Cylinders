% looping over different angles theta
r = 9/12/2;         % radius, ft
h = 51/12;          % height, ft
m = 113;            % mass, ft
theta = pi*[0.0557, 1/15, 1/12, 1/10, 1/8, 1/6, 1/5, 1/4];  % inclination
n = length(theta);
II_xbar = zeros(3,n);
for i = 1:n
    [I_xbar(:,i),II_xbar(:,i)] = get_config(r,h,theta(i));
end

% % horizontal distance between top of cylinders
% breadth = zeros(n,1);
% for i = 1:n
%     breadth(i) = II_xbar(1,i)+h/2*sin(theta(i))-(I_xbar(1,i)+h/2*sin(theta(i)));
% end

