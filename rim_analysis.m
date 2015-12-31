% rim_analysis.m
%-------------------------------------------------------------------------%
% This program will extract data from pairs of images using the MATLAB 
% ginput command

% Input/Required Data:
%   Camera Calibration Toolbox
%   Image directory
%   left images
%   right images

% Directory structure:
%     - image_directory
%         -calibration_dir
%         -data_dir


% Directories in Search Path:
% tool_box_calib - MATLAB Camera Calibration Toolbox
% lsge_matlab - set of Geometric Tools, used to fit a plane to data

% NOTE: Requires images be named according to '\trial_01_8862.tif' 

% Output:
% Cartesian coordinates 
% wavenumbers
%-------------------------------------------------------------------------%


% Initialization
%-------------------------------------------------------------------------%
clear 
% image_directory = 'D:\Drop Splash Data\11_01\';
tool_box_calib = 'C:\Users\roryhr\Documents\MATLAB\TOOLBOX_calib';
lsge_matlab = 'C:\Users\roryhr\Documents\MATLAB\lsge-matlab';

image_directory = 'C:\Users\roryhr\Documents\More Drop Splash Data - Processed\1_16\';
calibration_dir = 'Calibration\';
data_dir = 'data\';
image_data = 'image_data.mat';
Image_number = 1;
n_cell_width = 5;
hand = 'right'; % indicates which position the 8862 camera is in (left or right)

% Add toolboxes to search path
addpath(tool_box_calib, lsge_matlab)
%-------------------------------------------------------------------------%


% Load Calibration Data
% We are only loading in the variables we need 
stereo_calibration_results = 'Calib_Results_stereo.mat';
filename = strcat(image_directory, calibration_dir, stereo_calibration_results);
filevariables = {'om','T','fc_left','cc_left','kc_left', ...
    'kc_left','alpha_c_left','fc_right','cc_right','kc_right', ...
    'alpha_c_right', 'R', 'T', 'omc_left_1', 'Tc_left_1' };
load(filename, filevariables{:})


%% Image Acquisition 
%---------------------------------------------%
% hand indicates which position the 8862 camera is in (left or right)
[all_images] = image_acq(image_directory, data_dir , image_data, hand);

%% Data Acquisition 
%-------------------------------------------------------------------------%
% This loop allows the user to specify an image and pick corresponding pairs 
% of points in that image. Calibration results are used to determine the 
% location in cartesian space of each peak. 
% n = size(I_all_left);

while(Image_number ~= 0)
    
    % Choose image and frame to analyze
    %---------------------------------------------%
    dir(image_directory)
    valid_num = size(all_images.left,2);
    message = sprintf('Enter Image Number (1-%d, []=exit): ',valid_num);
    Image_number = input(message);
    
    if(isempty(Image_number)||(Image_number>valid_num))
        break
    end

    figure
    n = size((all_images.right(Image_number).image),2);
    m = ceil(n/n_cell_width);
    for i = 1:n
        subplot(m,n_cell_width,i) 
        subimage(all_images.left(Image_number).image(i).frame)
        set(gca,'xtick',[],'ytick',[])
        title(num2str(i))
    end
    
    message = sprintf('Enter frame number(1-%d, []=exit): ',n);
    frame_number = input(message);
    if(isempty(frame_number)||(frame_number>n))
        break
    end   
   
%%  	Pick points from left and right images
%-------------------------------------------------------------------------%
    [x_left, n] = pick_points(all_images, Image_number, frame_number, 'left');
    if(n==0)
        continue
    end
    message = sprintf('\nExtracted %d points from frame %d of: \n%s', ...
    	n, frame_number, all_images.left(Image_number).file_name);
    disp(message)

    [x_right, n] = pick_points(all_images, Image_number, frame_number, 'right');
   	if(n==0)
        continue
    end
    message = sprintf('Extracted %d points from frame %d of: \n%s', ...
    	n, frame_number, all_images.right(Image_number).file_name);
    disp(message)
%-------------------------------------------------------------------------%
    
    
%%     % Stereo Triangulation
%-------------------------------------------------------------------------%
  	[Xc_1_left,Xc_1_right] = stereo_triangulation(x_left',x_right',om,T, ...
 	fc_left, cc_left,kc_left,alpha_c_left,fc_right,cc_right, ...
  	kc_right,alpha_c_right); 
%% 


%%     % Plane fit
%-------------------------------------------------------------------------%
% Debug:
% Xc_1_left = [0,1,1; 0,-1,2; 0,1,9; 0,3,2]'

% Fit a plane in a least squares sense
% a - direction cosines of normal to best fit plane
[x0, a, d, normd] = lsplane(Xc_1_left');

% v = vector that gets rotated into a (by angle theta)
v = [0;0;1];

% u - vector in the fitted plane and the axis about which we rotate
u = cross(v,a);
u = u/norm(u);

theta = acos(dot(v,a));
% sin_theta = (1 - cos_theta^2 )^.5;

% v_rotation = v*cos_theta + cross(u,v)*sin_theta + z*(dot(u,v))*(1 - cos_theta);
[R] = rodrigues(-u.*theta);
XX = R * Xc_1_left;

%%
%     % Reference Frame transform
%     %---------------------------------------------%
% %     The first calibration image is our reference frame. We'll work in the
% %     left reference frame. 
%     % Rotation vector: omc_left_1 
%     % Translation vector:  Tc_left_1 
% %     These two variables are contained within Calib_Results_stereo.mat
%     % XXc = Rc_1 * XX + Tc_1 
% 
%     R = rodrigues(omc_left_1);
%     
%     Xc_1_t(1,:) = Xc_1_left(1,:)- Tc_left_1(1);
%     Xc_1_t(2,:) = Xc_1_left(2,:)- Tc_left_1(2);
%     Xc_1_t(3,:) = Xc_1_left(3,:)- Tc_left_1(3);
% 
%     XX = R\Xc_1_t;
%     figure
%     plot3(XX(1,:),XX(2,:),XX(3,:))
%     figure
%     scatter(XX(1,:),XX(2,:))
% 
%     
% %%     % Fit a circle in 3D and extract angles
% %-------------------------------------------------------------------------%
% guess_center =  mean(Xc_1_left')';
% guess_radius = mean(max(Xc_1_left') - min(Xc_1_left'))/2;
% [x0n, an, rn, d, e, f, sigmah, conv, Vx0n, Van, urn, ... 
%             GNlog, a, R0, R] = ls3dcircle(Xc_1_left', x0, a, guess_radius, 1, 1);
    
    
    

%%     % Fit a circle in 2D and extract angles
    %---------------------------------------------%
    % Taken from MATLAB example "Measuring the Radius of a Roll of Tape"
    % We begin with an equation of the form:
    % (x-xc)^2 + (y-yc)^2 = radius^2,  where (xc,yc) is the center. 
    % x^2 - 2*x*xc + xc^2 + y^2 - 2*y*yc + yc^2 - radius^2 = 0
    % In terms of parameters a,b,c:
    % x^2 + y^2 + a*x + b*y + c = 0
    % where a = -2*xc, b = -2*yc, and c = xc^2 + yc^2 - radius^2
    
%     We assume that XX(3,:) = const, and project points onto x-y plane
%     The backslash operator solves for a,b,c in a least squares sense
    abc = [XX(1,:)' XX(2,:)' ones(length(XX(1,:)),1)] \ -(XX(1,:)'.^2+XX(2,:)'.^2);
    a = abc(1); b = abc(2); c = abc(3);

    % calculate the location of the center and the radius
    xc = -a/2;
    yc = -b/2;
    radius  =  sqrt((xc^2+yc^2)-c);

    % Move XX to the origin
    XX2(1,:)= XX(1,:)' - xc;
    XX2(2,:)= XX(2,:)' - yc;

    % Calculate angles using atan2
    theta = atan2(XX2(2,:),XX2(1,:));
    % figure, scatter(XX(1,:),XX(1,:))

%%     % Wavenumber extraction
%-------------------------------------------------------------------------%

% Complex Order Parameter calculation
%-------------------------------------------------------------------------%
    N = length(theta);
    r = zeros(N,1);

    for m = 1:N;  
        r(m) = sum( exp(m*1i * theta))/N;
    end
    
    m = 1:N;
    figure, scatter(m,abs(r))
    xlabel('m')
    ylabel('r(m) - order amplitude')
%-------------------------------------------------------------------------%

[C,I] = max(r);
k_test = m(I); % Take the largest r value as the most likely wavenumber

[k, points] = wavenumber3(theta', k_test);


%%     % Save results
    %---------------------------------------------%
    % Here we store the data points extracted, the fitted circle, and the
    % extracted angles
    if(~isdir(strcat(image_directory, data_dir)))
        mkdir(strcat(image_directory, data_dir))
    end
    filename = sprintf('%strial_%d_frame_%d',strcat(image_directory, data_dir),Image_number,frame_number);   
    save(filename)
    clear XX2

%     fid = fopen(strcat(image_directory, data_dir, filename),'at+');
%     fprintf(fid, '\n Image Number: %d \t Frame Number: %d \t \n r=%f, xc=%f, yc=%f ', ...
%         Image_number, frame_number, radius, xc, yc);
%     fprintf(fid, '\n X \t\t Y \t\t Z \t\t Theta');
%     fprintf(fid, '\n %f \t %f \t %f \t %f', XX2(1,:)', XX2(2,:)', XX(3,:)', theta);
%     fprintf(fid, '\n -------------------------------------------------------- \n');
%     fclose(fid);
%     
end
