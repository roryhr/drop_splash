function [] = complex_order_parameter()
%-------------------------------------------------------------------------%
%% DEBUG
k1 = 9;
k2 = 14;
phi_range = 1.2;
error = 0.1;
[theta] = gen_pattern(k1, k2, phi_range, error);
%-------------------------------------------------------------------------%
% figure
% scatter(cos(theta),sin(theta))

function [theta] = gen_pattern(k1, k2, phi_range, error)
phi = phi_range*rand(1); % Mean phase shift of 45 degrees
a = linspace(-pi, pi, k1+1) + phi ;
b = linspace(-pi, pi, k2+1);

theta = unique([a b])';
theta = theta + error*rand(size(theta))- error/2;  %Add zero mean error

theta(theta > pi) = theta(theta > pi) - 2*pi;
theta(theta < -pi) = theta(theta < -pi) + 2*pi;
theta = sort(theta);
theta = theta(abs([diff(theta)' (theta(end)+ theta(1))] ) >= error/2); % Consoladate near points
% figure, scatter(cos(theta), sin(theta))
end

N = length(theta);
r = zeros(N,1);

for m = 1:N;
    r(m) = sum( exp(m*1i * theta))/N;
end
    
figure
m = 1:N;
scatter(m,abs(r))
xlabel('m')
ylabel('r(m) - order amplitude')
end
