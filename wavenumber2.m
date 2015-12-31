function [k, points] = wavenumber2(theta)
% wavenumber.m
% [k, points] = wavenumber(theta)


% Input:
% theta: an (n,1) vector consisting of angles between (-pi, pi)

% Output:
% Wavenumber - vector of wavenumbers
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
% n_iters = 2; % Number of iterations, ie expected number of wavenumbers
d2 = 0;
min_length = 3; % Min number of points required 
%-------------------------------------------------------------------------%

    
% First form combinations of indicies of theta 
c = combnk(1:length(theta),2); 
% Then evaluate and compute differences  
% By construction, these differences of theta are in the range of 0-2*pi 
d1 = diff(theta(c),1,2);


% Make a hist plot of finite differences
%-------------------------------------------------------------------------%
d1(d1==0)= .0001; % Remove sigularities
k1 = 2:floor(1.5*length(theta));

n = mod_histogram(d1,k1, tol3);
% figure1 = figure;
% axes1 = axes('Parent',figure1,'YGrid','on','XMinorTick','on');
% box(axes1,'on');
% hold(axes1,'all');
% ind2 = mod(ind,1) == 0 & ind > k_lower_thresh ;
% bar(k1,n)
%     axes('XTick', ind)
% We will now extract points and verify
[C I] = max(n);
%-------------------------------------------------------------------------%


[x_n1,  theta_1] = wave_fit(theta, d1, c, k1(I), tol3);

theta_diff = setdiff(theta, theta_1);

if length(theta_diff) > min_length
    % First form combinations of indicies of theta 
    c2 = combnk(1:length(theta_diff),2); 

    d2 = diff(theta(c2),1,2);
    k2 = 2:floor(1.5*length(theta_diff));
    n2 = mod_histogram(d2,k2, tol3);
    n2(k2==k1(I)) = 0 ; % Remove k1(I) from our search
    n2(k2==k1(I)/2) = 0 ; % 1st harmonic
    
    figure1 = figure;
    axes1 = axes('Parent',figure1,'YGrid','on','XMinorTick','on');
    box(axes1,'on');
    hold(axes1,'all');
    bar(k2,n2)
    
    
    
    [C2 I2] = max(n2);
    [x_n2,  theta_2] = wave_fit(theta_diff, d2, c2, k2(I2), tol3);

    figure
    scatter(x_n1, theta_1), hold on
    [k1 u1 p1] = linear_k_fit(x_n1, theta_1);
    plot(x_n1, polyval(p1, x_n1))

    scatter(x_n2, theta_2)
    [k2 u2 p2] = linear_k_fit(x_n2, theta_2);
    plot(x_n2, polyval(p2, x_n2)) 
    
%-------------------------------------------------------------------------% 
message = sprintf('\nk''s: %.1f +- %.1f, %.1f +- %.1f, ', k1, u1, k2, u2);            
disp(message)
% message = sprintf('%2.2f  ', theta_prov3(ind));            
% disp(message)
end


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
        elseif length(ind2) < k_prov/2
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
    
%     figure, scatter( x_n2,  theta_prov2 )


% Given our line defined by k3, we'll find points in theta which are
% closest to the line
%-------------------------------------------------------------------------%
%     theta_prov3(:,1) = 1:k_prov;
%     for i = 1:length(theta_prov3)                 
%         [res3(i) I3(i)] = min(abs(polyval(p3,theta_prov3(i,1)) - theta));                  
%     end 
%     theta_prov3(:,2) =  theta(I3);

% % Remove any "jumps" in our plot of theta      
%     if sum(diff(theta_prov3(:,2)) < 0.5*dth)            
%         a3 = round([diff(theta_prov3)', abs(theta(1)+theta(length(theta)))]'/dth)-1;
%         x_n = x_n + cumsum(a2) - a2;
%     end
%     
%     figure, scatter(theta_prov3(:,1),theta_prov3(:,2))
%     [k4 u4] = linear_k_fit(x_n2(ind), theta_prov3(ind));

   
end

function [n] = mod_histogram(theta, k, epsilon)
% theta(theta==0) = 0.00001; %Remove possible singular points 
n = zeros(length(k), 1);
for i=1:length(k)
    dth = 2*pi/k(i);
    n(i) = sum( abs( mod(theta, dth)./dth - 1 ).* ceil(theta./dth) < epsilon );
end 
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
sig_y = (sum( (y - polyval(p, x)).^2 )/(N-2))^.5;
del = N * sum(x.^2) - (sum(x))^2;
sig_a = sig_y * ((sum(x.^2))/del)^.5;

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