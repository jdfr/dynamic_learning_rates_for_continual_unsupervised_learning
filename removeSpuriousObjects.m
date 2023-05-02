function bw_out = removeSpuriousObjects(bw, min_area)
% bw_out = removeSpuriousObjects(bw, min_area)
% The objects in the image 'bw' with few pixels (less than 'min_area') will
% be removed in the output image. 

% R.M.Luque and Ezequiel Lopez-Rubio -- June 2011

[L,num_blobs] = bwlabel(bw);
        
bw_out = zeros(size(bw,1),size(bw,2));
for i=1:num_blobs 
    Ob = (L == i);
    area = bwarea(Ob);
    if area > min_area
        bw_out = bw_out + Ob;
    end
end
