function [X_p_out] = point_line(X_p_in, fc_left, fc_right, cc_left, ...
    cc_right, kc_left, kc_right, alpha_c_left, alpha_c_right, om, T, hand)
% point_line.m
% [X_p2] = point_line([50;50],fc_left,cc_left,kc_left,alpha_c_left, om, T, 'left')
    
%INPUT: 
%       X_p_in: (2xN) vector of (x,y) pixel coordinates
%       (om,T): Rigid motion parameters between world coordinate frame and camera reference frame
%       om: rotation vector of right wrt left(3x1 vector) 
%       T: translation vector to right wrt left (3x1 vector)
%       fc: camera focal length in units of horizontal and vertical pixel units (2x1 vector)
%       cc: principal point location in pixel units (2x1 vector)
%       k: Distortion coefficients (4x1 vector)
%       alpha: Skew coefficient between x and y pixel

%OUTPUT: 
% X_p_out = (x,y) pixel coordinates in the other camera

    Z_c = 100:.5:1000;      % define Z_c here
%     R = rodrigues(om); % generate the corresponding rotation matrix (the identity matrix) 
    
    if strcmpi( hand, 'left') 
%         We are given points in the left image and wish to correspond
%         these to lines in the right image
        X_n_left = normalize(X_p_in',fc_left,cc_left,kc_left,alpha_c_left);
        
        XX_c_left = X_n_left * Z_c;
        XX_c_left(3,:) = Z_c;
        
%         XX_c_right = R * XX_c_left +  repmat(T,[1 size(XX_c_left,2)]);
        
        [X_p_right,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk,dxpdalpha] = ...
            project_points2(XX_c_left, om, T, fc_right,cc_right,kc_right,alpha_c_right);
        X_p_right = X_p_right';
        
        [row col] = find( X_p_right(:,:) > 0 & X_p_right(:,:) < 769 ); % Restrict to pixel area
        X_p_out = X_p_right(row,:);
        
    elseif strcmpi( hand, 'right')      
        return
    else
        disp('Must select either left or right hand')
        return
    end
        


end