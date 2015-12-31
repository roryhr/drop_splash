function [k, points] = wavenumber()
% wavenumber.m
% [k, points] = wavenumber(theta)


% Input:
% theta: an (n,1) vector consisting of angles between (-pi, pi)

% Output:
% Wavenumber - vector of wavenumbers
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%% DEBUG
clear
% theta =[0; 0.4487; 1.5708; 2.5462; 3.1416 ;-1.6442; -1.5708];
a = linspace(-pi, pi, 10) + 2*rand(1) ;
b = linspace(-pi, pi, 15);

theta = unique([a b])';
theta = theta + 0.01*rand(size(theta))-0.005;
k_lower_thresh = 5; % Lower threshold for wavenumber
theta(theta > pi) = theta(theta > pi) - 2*pi;
theta(theta < -pi) = theta(theta < -pi) + 2*pi;
theta = sort(theta);
theta = theta(abs([diff(theta)' (theta(length(theta))+ theta(1))] ) > .06); % Consoladate near points
figure, scatter(cos(theta), sin(theta))
%-------------------------------------------------------------------------%


%% Specify tol, which is the histogram bin width
%-------------------------------------------------------------------------%
tol = .5; % Bin width
tol2 = 0.09; % Rejection of points which are too close
tol3 = 0.08; % Find points corresponding to a wavenumber
n_iters = 2; % Number of iterations, ie expected number of wavenumbers
d2 = 0;
%-------------------------------------------------------------------------%

    
% First form combinations of indicies of theta 
c = combnk(1:length(theta),2); 
% Then evaluate and compute differences  
d1 = diff(theta(c),1,2);
 
%     Put differences of theta into the range of 0-pi 
%-------------------------------------------------------------------------% 
d2 = mod2pi(d1);
%-------------------------------------------------------------------------%


%     First make a hist plot of finite differences
%-------------------------------------------------------------------------%
d2(d2==0)= .001; % Remove sigularities
k1 = 2*pi./d2;
[n ind]  = hist(k1, 1:tol:(length(theta)+1));
figure1 = figure;
axes1 = axes('Parent',figure1,'YGrid','on','XMinorTick','on');
box(axes1,'on');
hold(axes1,'all');
ind2 = mod(ind,1) == 0 & ind > k_lower_thresh ;
bar(ind(ind2),n(ind2))
%     axes('XTick', ind)
% We will now extract points and verify
  
ind2 = find( n > mean(n(n>0))/2 );
k_prov = ind(ind2);
k_prov = k_prov(mod(k_prov,1) < tol/2 & k_prov > k_lower_thresh);

for l = 1:length(k_prov)
    wave_fit(theta, d2, k_prov(l))
end


function [] = wave_fit(theta, d2, k_prov)
    %% Find appropriate angles
    ind = abs(d2 - 2*pi/k_prov) < tol3 ; %Find these wavenumbers
    ind2 = unique(c(ind,:)); %These are the indicies of the original theta vector

    % Generate a provisional theta vector, sorted in ascending order
    theta_prov = sort(theta(ind2));
    x_n = (1:length(theta_prov))'; %Basically, the indicies of theta_prov
    [k2 u2] = linear_k_fit(x_n, theta_prov);


%-------------------------------------------------------------------------%
    % Evaluate another difference of theta_prov to ensure that we aren't
    % duplicating a theta value which would throw off the index
    theta_diff_prov = abs(diff(theta_prov));
    theta_prov2 = theta_prov;
  
    for i = 1:length(theta_diff_prov)
        if theta_diff_prov(i) < tol2
            message = sprintf('Points too close, merging values: %02.2f and %02.2f of theta_prov', ...
                theta_prov(i), theta_prov(i+1));
            disp(message)
%             theta_prov2(i+1) = [];
            theta_prov2(i+1,1) = (theta_prov(i) + theta_prov(i+1))/2; 
            index_flag(i) = 1;
        else
            index_flag(i) = 0;           
        end
    end
    theta_prov2(index_flag==1) = [];

    
%     Compute k a second time
%-------------------------------------------------------------------------%
    x_n = (1:length(theta_prov2))';
%     figure
%     scatter(x_n,theta_prov2)
    [k3 u3 p] = linear_k_fit(x_n, theta_prov2);
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
    %% Recompute k a third time
%     This time we'll associate the points in theta with that closest
%     to our fitted line. Then recalculate.
%     theta_prov3 = theta_prov2;
    x_n2 = (1:length(theta_prov2))';
    for i = 1:length(x_n2)         
        [C I] = min(abs(polyval(p,x_n2(i)) - theta ));
        res(i,1) = polyval(p,x_n2(i)) - theta(I);
        theta_prov3(i,1) =  theta(I);
    end
    
    ind = res < mean(res) + std(res);

%-------------------------------------------------------------------------%
    % Fit a first order polynomial
%     figure
%     p = polyfit( x_n2(ind), theta_prov3(ind), 1);
    figure, scatter( x_n2(ind),  theta_prov3(ind) )
    [k4 u4] = linear_k_fit(x_n2(ind), theta_prov3(ind));

    
%-------------------------------------------------------------------------% 
message = sprintf('\nInitial guess %.1f, iterated k''s: %.1f +- %.1f, %.1f +- %.1f, %.1f +- %.1f\nAssociated theta values: ', ...
    k_prov, k2,u2, k3, u3, k4, u4);            
disp(message)
message = sprintf('%2.2f  ', theta_prov3(ind));            
disp(message)

%-------------------------------------------------------------------------%
%     %% Remove the points contained in theta_prov2 from theta
    theta_diff = setdiff(theta,theta_prov3);
    clear ind ind2 ind3 theta_prov theta_prov2 theta_prov3 d1 x_n2 x_n index_flag res
%     theta = theta_diff;
%     %iterate...
% %-------------------------------------------------------------------------%
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