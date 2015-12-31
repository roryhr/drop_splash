function [x_points, n] = pick_points(all_images, Image_number, frame_number, hand)
% Function that collects points from an image via curser input. 

    figure
    axes
    if strcmpi( hand, 'left') 
%         im_mean = mean(mean(all_images.left(Image_number).image(frame_number).frame));
%         Use im_mean to apply a window to enhance contrast of left image
%         imshow(all_images.left(Image_number).image(frame_number).frame, [im_mean - 30, im_mean+30])
        imshow(all_images.left(Image_number).image(frame_number).frame)
    elseif strcmpi( hand, 'right')
%         [X_p_out] = point_line(x_left(1,:), fc_left, fc_right, cc_left, ...
%     cc_right, kc_left, kc_right, alpha_c_left, alpha_c_right, om, T, hand);

        imshow(all_images.right(Image_number).image(frame_number).frame) 
%         hold on, scatter(X_p_out(:,1), X_p_out(:,2), '.y')
    else
        disp('Must select either left or right image from which to pick points')
        return
    end
        
    hold on
    % Initially, the list of points is empty.
    x_points = [];
    n = 0;
    % Loop, picking up the points.
    % Left mouse button picks points
    % Right mouse button picks last point
    but = 1;
    while but == 1
        [xi,yi,but] = ginput(1);
        if but==1
            n = n+1;
            plot(xi,yi,'ro','MarkerSize',5)
            text(xi+5,yi+9, num2str(n, '%.0f'),'FontSize',12,'Color','y')
            x_points(n,:) = [xi,yi];
        end        
    end
end