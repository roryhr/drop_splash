function [x_n1, theta_1, k1, u1, p1] = wavenumber3(theta, k_test)
% wavenumber.m
% [k, points] = wavenumber(theta)


% Input:
% theta: (n,1) vector consisting of angles between (-pi, pi)
% k_test: (n,1) vector consisting of wavenumbers to analyze

% Output:
% x_n1: "x-coordinates"
% theta_1: fitted angles 
% k1: wavenumber calculated from slope of fitted line
% u1: uncertainty in k1
% p1: fitted polynomial
%-------------------------------------------------------------------------%


% %-------------------------------------------------------------------------%
% %% DEBUG
% k1 = 8;
% k2 = 11;
% phi_range = 1.2;
% error = 0.1;
% [theta] = gen_pattern(k1, k2, phi_range, error);
% %-------------------------------------------------------------------------%


%% Specify tol, which is the histogram bin width
%-------------------------------------------------------------------------%
tol = .5; % Bin width
tol2 = 0.09; % Rejection of points which are too close
tol3 = 0.08; % Find points corresponding to a wavenumber

d2 = 0;
min_length = 3; % Min number of points required 
%-------------------------------------------------------------------------%

    
% First form combinations of indicies of theta 
% c = combnk(1:length(theta),2); 
c = combnk2(1:length(theta)); 
% Then evaluate and compute differences  
% By construction, these differences of theta are in the range of 0-2*pi 
d1 = diff(theta(c),1,2);


% Make a hist plot of finite differences
%-------------------------------------------------------------------------%
d1(d1==0)= 0.0001; % Remove sigularities
% k1 = 2:floor(1.5*length(theta));

[x_n1,  theta_1] = wave_fit(theta, d1, c, k_test(1), tol3);

figure
scatter(x_n1, theta_1), hold on
    [k1 u1 p1] = linear_k_fit(x_n1, theta_1);
    plot(x_n1, polyval(p1, x_n1))

%-------------------------------------------------------------------------% 
message = sprintf('\nk: %.1f +- %.1f ', k1, u1);            
disp(message)


end

function [x_n2,  theta_prov2] = wave_fit(theta, d2, c, k_prov, epsilon)
% Find appropriate angles
%-------------------------------------------------------------------------%
    dth = 2*pi/k_prov;
    x=0;
    for i = 1:10 %Apply a fixed number of iterations to our "adaptive filter"
        ind = abs(mod(d2, dth)./dth - 1 ).* abs(ceil(d2./dth)) < epsilon ; %Find these wavenumbers
        ind2 = unique(c(ind,:)); %These are the indicies of the original theta vector
        
        if length(ind2) > k_prov
            % We have too many angles, tighten epsilon
            x = x + 1;
        elseif length(ind2) <= (k_prov-1)
            % We have less than half the expected angles, loosen epsilon
            x = x - 1;
        else
            break
        end
        epsilon = epsilon *(.9)^x ;
    end
%-------------------------------------------------------------------------%
    
% Generate a provisional theta vector, sorted in ascending order
%-------------------------------------------------------------------------%
    theta_prov = sort(theta(ind2));
    x_n = (1:length(theta_prov))'; %Basically, the indicies of theta_prov    
    % Remove any "jumps" in our plot of theta   
%     figure, scatter(x_n, theta_prov)
    if sum(diff(theta_prov) > 1.5*dth)            
        a2 = round([diff(theta_prov)', abs(theta(1)+theta(length(theta)))]'/dth)-1;
        x_n = x_n + cumsum(a2) - a2;
    end
        
%     figure, scatter(x_n, theta_prov)
    [k2 u2 p2] = linear_k_fit(x_n, theta_prov);

    
% Remove duplicated angles in theta_prov by selecting the one closest to
% our fitted line
%-------------------------------------------------------------------------%
I1 = find(diff(x_n) == 0);
if ~isempty(I1)
    theta_prov2 =  theta_prov;
        for i = 1:length(I1)         
            [res(i) I(i)] = min(abs(polyval(p2,x_n(I1(i))) - theta_prov ));           
        end 
    theta_prov2(I1) =  theta_prov(I);
    theta_prov2(I1+1) = [];
    x_n2 = x_n;
    x_n2(I1) = [];
%-------------------------------------------------------------------------%


% Compute k a second time, noting that we have removed duplicated points in
% theta_prov
%-------------------------------------------------------------------------%
    [k3 u3 p3] = linear_k_fit(x_n2, theta_prov2);
%-------------------------------------------------------------------------%
else
%     We haven't made any of the previous calculations, so simply move on
%     by letting var2 = var
    k3 = k2; u3 = u2; p3 = p2; x_n2 = x_n;
    theta_prov2 = theta_prov;
end 
end


function [c] = combnk2(x)
% Duplicates functionality of Matlab function combnk(x,2)
index = 0;
N = length(x);
c = zeros(factorial(N)/(factorial(N-2)* 2),2);
for i = 1:(N-1)
    for j = (i+1):N
        index = index + 1;
        c(index,1) = i;
        c(index,2) = j;
    end 
end
c = flipud(c);
end

function [angle_out] = mod2pi(angle_in)
    for i = 1:length(angle_in)
        if angle_in(i)< -pi 
            angle_out(i) = angle_in(i) + 2*pi;
        elseif angle_in(i)> pi
            angle_out(i) = angle_in(i) - 2*pi; 
        else
            angle_out(i) = angle_in(i);    
        end            
        if abs(angle_out(i))<=pi       
            angle_out(i) = abs(angle_out(i));           
        end       
    end
end

function [k uncertainty p] = linear_k_fit(x, y)
% Fit a first order polynomial
p = polyfit( x, y, 1);
k = 2*pi/p(1);
% Uncertainty Analysis
N = length(y);
sig_y = sqrt(sum( (y - polyval(p, x)).^2 )/(N-2));
del = N * sum(x.^2) - (sum(x))^2;
sig_a = sig_y * sqrt( (sum(x.^2))/del );

uncertainty = 2*pi*sig_a;
end

function [theta] = gen_pattern(k1, k2, phi_range, error)
phi = phi_range*rand(1); % Mean phase shift of 45 degrees
a = linspace(-pi, pi, k1+1) + phi ;
b = linspace(-pi, pi, k2+1);

theta = unique([a b])';
theta = theta + error*rand(size(theta))- error/2;  %Add zero mean error
% k_lower_thresh = 5; % Lower threshold for wavenumber
theta(theta > pi) = theta(theta > pi) - 2*pi;
theta(theta < -pi) = theta(theta < -pi) + 2*pi;
theta = sort(theta);
theta = theta(abs([diff(theta)' (theta(length(theta))+ theta(1))] ) ~= 0); % Consoladate near points
% figure, scatter(cos(theta), sin(theta))
end