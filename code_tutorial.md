# CNC-SMART for TEM image analysis
## To select and upload TEM images
```
clear all,clc
[path,user_cance] = imgetfile('MultiSelect','on');
nfiles = length(path);
ims = cell(nfiles,1);
for i = 1:nfiles
    ims{i} = imread(path{i});
end
%%
```
## To enter pixel length information (using command window)
```
pixel_length = input("What is the original value? ");
```
## Initiation of variables (no user input needed)
```
allL = [];
allW = [];
allL_max = [];
allW_max = [];
border = cell(nfiles,1);
isolated = cell(nfiles,1);
overlaps = cell(nfiles,1);
isolatednew = cell(nfiles,1);
parallels = cell(nfiles,1);
count = 0;
filter1 = [1;-1];
filter2 = [1,-1];
for i = 1:nfiles
 	img = (ims{i});
	img = colfilt(img,[3 3],'sliding',@max);
  med_size = 1+2*floor(0.125*15*1/pixel_length);
  img = imadjust(img);
	lsharp = 0.5;
  img1 = imsharpen(img,'radius',lsharp);
	img1 = medfilt2(img1,[med_size med_size],'symmetric');
	nb_size = 1+8*floor(15*0.5/pixel_length);
	T = adaptthresh(img1,0.5,'ForegroundPolarity','bright','NeighborhoodSize',[nb_size nb_size]);
	bwi = imbinarize(img1,T);
	bwi1 = imclearborder(bwi);
	bwb = bwi-bwi1;
	bwi = bwareaopen(bwi1,floor(10*5/pixel_length));
	bwi = imerode(bwi,strel('disk',2));
	bwi = imdilate(bwi,strel('disk',1));
	bwi = bwareaopen(bwi,floor(10*5/pixel_length));
	bw2 = bwareaopen(bwi,20*floor(20*1/pixel_length));
	bwall = bw2+bwb;
	idbwi = find(bwall==0);
	idcnc = find(bwall == 1);
	imgback = (ims{i});
	noise_std(i) = std(double(imgback(idbwi)));
	contrast(i) = mean(double(imgback(idbwi)))-mean(double(imgback(idcnc)));
	...
	
end
```
