function [all_images] = image_acq(image_directory, data_dir , image_data, hand)
% Returns a struct all_images containing all the left and right images 
% hand indicates which position the 8862 camera is in (left or right)

% NOTE: Requires images be named according to '\trial_01_8862.tif' - the
% trial number is read from elements [end-10:end-9]

% Image Acquisition 
%---------------------------------------------%
filename = strcat(image_directory, data_dir , image_data);
message = sprintf('\nImage directory: %s', image_directory);
disp(message)
if(isempty(input('Load previously stored images([]=yes): ', 's')))
    load(filename)
    return
else  
    if strcmpi( hand, 'left')
        image_list_left = dir(strcat(image_directory,'*8862*.tif'));
        image_list_right = dir(strcat(image_directory,'*6444*.tif'));
    elseif strcmpi( hand, 'right')
        image_list_left = dir(strcat(image_directory,'*6444*.tif'));
        image_list_right = dir(strcat(image_directory,'*8862*.tif'));
    end
    
   
    % Check the directory for all .tif files and read in a stack for each .tif file
    % Loop through all the files in the folder

    for j = 1:length(image_list_left)
        filename_left = strcat(image_directory, image_list_left(j).name);
        filename_right = strcat(image_directory, image_list_right(j).name);
        image_info_left = imfinfo(filename_left);
        image_info_right = imfinfo(filename_right);
        
        trial_number = str2num(filename_left(end-10:end-9));
        
        all_images.left(trial_number).file_name = filename_left;
        all_images.right(trial_number).file_name = filename_right;
        
        for i = 1:size(image_info_left,1)
%             I_im_left(:,:,i) = imread(filename_left,i);
%             I_im_right(:,:,i) = imread(filename_right,i);
            all_images.left(trial_number).image(i).frame = imread(filename_left,i);
            all_images.right(trial_number).image(i).frame = imread(filename_right,i);
        end

    end 
    
end

if(isdir(strcat(image_directory, data_dir)))
    save(filename, 'all_images')  
else
    mkdir(strcat(image_directory, data_dir))
    save(filename, 'all_images')  
end

end
