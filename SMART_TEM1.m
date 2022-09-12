%% to select and upload the TEM images
%clear all,close all, clc
[path,user_cance] = imgetfile('MultiSelect','on');
nfiles = length(path);
ims = cell(nfiles,1);
for i = 1:nfiles
    ims{i} = imread(path{i});
end
%% the first input to be entered by the user
input1 = inputdlg('Enter the pixel length (nm):');
pixel_length = str2num(input1{1});
