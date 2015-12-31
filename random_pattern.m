% random_pattern.m
n = 111; %Number of randomly generated points 
theta = 2*pi*rand(n,1);
theta = sort(theta);

lambda = 2*pi/n;

theta2 = theta(1) + [0:n-1]'*lambda;
theta2= mod(theta2, 2*pi);




figure, scatter(cos(theta), sin(theta))
hold on
scatter(cos(theta2), sin(theta2),'x')
figure
hist((theta - theta2)/lambda)