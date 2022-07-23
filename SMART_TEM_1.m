%% the first parameter to be entered by the user:
pixel_length= 0.33;
%% to select and upload the TEM images
[path,user_cance] = imgetfile('MultiSelect','on');
nfiles = length(path);
ims = cell(nfiles,1);
for i = 1:nfiles
    ims{i} = imread(path{i});
end