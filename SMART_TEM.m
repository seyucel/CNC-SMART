%% to select and upload the TEM images
%clear all,close all, clc
[path,user_cance] = imgetfile('MultiSelect','on');
nfiles = length(path);
ims = cell(nfiles,1);
for i = 1:nfiles
    ims{i} = imread(path{i});
end
%% the first input to be entered by the user:
input1 = inputdlg('Enter the pixel length (nm):');
pixel_length = str2num(input1{1});
%% initial parameters // can be adjusted later
lsmooth = 5;
lsharp = 1;
%
allL = [];
allW = [];
allL_max = [];
allW_max = [];
border = cell(nfiles,1);
pre_bw = cell(nfiles,1);
isolated = cell(nfiles,1);
overlaps = cell(nfiles,1);
isolatednew = cell(nfiles,1);
parallels = cell(nfiles,1);
filter1 = [1;-1];
filter2 = [1,-1];
nb_size = 1+2*floor(30/pixel_length);
med_size = 1+2*lsmooth; %+ 2*floor(0.075*15*1/pixel_length);
% segmentation/image quality check for the first image
i = 1;
img = ims{i};
ids = find(img==0);
id_scale_bar = find(img == 255);
img(id_scale_bar) = round(mean(img(:)));
img = imadjust(img);
img1 = imsharpen(img,'radius',lsharp);
img1 = medfilt2(img1,[med_size med_size],'symmetric');
T = adaptthresh(img1,0.5,'ForegroundPolarity','bright','NeighborhoodSize',[nb_size nb_size]);
bwi = imbinarize(img1,T);
bwi1 = imclearborder(bwi);
bwb = bwi-bwi1;
bwi = bwareaopen(bwi1,floor((1/pixel_length)^2)*150);%50/pixel_length
bwi = imerode(bwi,strel('disk',2));
bwi = imdilate(bwi,strel('disk',2));
bwi = bwareaopen(bwi,floor((1/pixel_length)^2)*150,4);
bwi = bwareaopen(bwi,floor((1/pixel_length)^2)*150,4);
bwi = imerode(bwi,strel('disk',2));
bwi = imdilate(bwi,strel('disk',2));
bwi = bwareaopen(bwi,floor((1/pixel_length)^2)*150,4);

s = regionprops(bwi, 'all');
c = bwconncomp(bwi,8);
mja = [s.MajorAxisLength];
mna = [s.MinorAxisLength];
ar = mja./mna;
bw2 = ismember(labelmatrix(c), find(mja>(22/pixel_length)));

bwb = bwareaopen(bwb,floor((1/pixel_length)^2)*100);
bwall = bw2+bwb;
bwall(1925:end,1:222) = 0;
idbwi = find(bwall==0);
idcnc = find(bwall == 1);
imgback = (ims{i});
border{i} = bwb;
pre_bw{i} = bw2;

noise_std(i) = std(double(imgback(idbwi)));
contrast(i) = mean(double(imgback(idbwi)))-mean(double(imgback(idcnc)));
cnc_percentage(i) = length(idcnc)/numel(bwall);

if mean(noise_std) > 0.15*range(double(imgback(:))) && mean(contrast) < 0.075*range(double(imgback(:)))
    fprintf('WARNING: noise or contrast might affect the particle identification and measurements; \n \n try to increase med_size by 2 unit or decrease the l_sharp by 0.5 unit \n \n')
end

figure,subplot(121),imshow(imgback),subplot(122),imshow(bwall)
mkdir results\image_analysis_figures
ans1 = questdlg('Apply the same segmentation for the rest of the images','Segmentation quality check:','Yes','No','Yes');
if length(ans1) < 3
    fprintf('Change the lsharp and/or lsmooth values (lines 12-13) \n');
    return
else
    for i = 1:nfiles
        img = (ims{i});
        ids = find(img==0);
        id_scale_bar = find(img == 255);      
        img(id_scale_bar) = round(mean(img(:)));
        img = imadjust(img);
        img1 = imsharpen(img,'radius',lsharp);  
        img1 = medfilt2(img1,[med_size med_size],'symmetric');
        T = adaptthresh(img1,0.5,'ForegroundPolarity','bright','NeighborhoodSize',[nb_size nb_size]);
        
        bwi = imbinarize(img1,T);
        %bwi = ~bwi;
        bwi1 = imclearborder(bwi);
        bwb = bwi-bwi1;      
        bwi = bwareaopen(bwi1,floor((1/pixel_length)^2)*150);%50/pixel_length
        bwi = imerode(bwi,strel('disk',2));       
        bwi = imdilate(bwi,strel('disk',2));
        bwi = bwareaopen(bwi,floor((1/pixel_length)^2)*150,4);  
        bwi = bwareaopen(bwi,floor((1/pixel_length)^2)*150,4);
        bwi = imerode(bwi,strel('disk',2));       
        bwi = imdilate(bwi,strel('disk',2));
        bwi = bwareaopen(bwi,floor((1/pixel_length)^2)*150,4);
        
        s = regionprops(bwi, 'all');
        c = bwconncomp(bwi,8);
        mja = [s.MajorAxisLength];
        mna = [s.MinorAxisLength];
        ar = mja./mna;
        bw2 = ismember(labelmatrix(c), find(mja>(22/pixel_length)));
        
        bwb = bwareaopen(bwb,floor((1/pixel_length)^2)*100);
        bwall = bw2+bwb;
        bwall(1925:end,1:222) = 0;
        idbwi = find(bwall==0);
        idcnc = find(bwall == 1);
        imgback = (ims{i});
        border{i} = bwb;
        pre_bw{i} = bw2;     
        
        noise_std(i) = std(double(imgback(idbwi)));
        contrast(i) = mean(double(imgback(idbwi)))-mean(double(imgback(idcnc)));
        cnc_percentage(i) = length(idcnc)/numel(bwall);
        
        f = figure('visible','off');
        subplot(121),imshow(imgback,[]),subplot(122),imshow(bwall)
        pname = path(i);
      
        name_ex = ['results\image_analysis_figures\', 'img',num2str(i),'_segmented_duo''.png'];
     saveas(f,name_ex)
        

    end
end
avg_noise = mean(noise_std)
avg_contrast = -mean(contrast)
avg_CNC_area_fraction = mean(cnc_percentage)
%% Grouping step; user can adjust the selection criteria:
max_minor = 15;
min_major = 15;
max_major = 300;
ar_min = 2.5;
% grouping check for the first image
i = 1;
k=0;
    bw2 = pre_bw{i}; 
    s = regionprops(bw2, 'all');
    c = bwconncomp(bw2,8);
    for oi = 1:length(s)
        convimmini = s(oi).ConvexImage;
        convimmini1 = imrotate(convimmini,90-s(oi).Orientation);
        smini = regionprops(convimmini1,'BoundingBox');
        bboxvalues = smini.BoundingBox;
        mna(oi) = bboxvalues(3);
        mja(oi) = bboxvalues(4);
    end
    ar = mja./mna;
    bw3 = ismember(labelmatrix(c), find(ar>ar_min & mna < (max_minor/pixel_length) & mja > (min_major/pixel_length) & mja < (max_major/pixel_length)));
    bwoverlap1 = logical(bw2-bw3);
    bw3 = bwmorph(bw3,'thin',1);
    bwoverlap1 = bwareaopen(bwoverlap1,floor(numel(bwi)*0.0004));
    isolated{i} = bw3;
    overlaps{i} = (bwoverlap1);
    clusters = bwoverlap1;
    ss = regionprops(bw3,'all');
    cc = bwconncomp(bw3,8);
    image = ims{i};
    coef1 = [];
    coefnew = [];
    coefnewnew = [];
    coefnew_a = [];
    coefnew_b = [];
        
    for j = 1:length(ss)
        angle = ss(j).Orientation;
       obj_img = ismember(labelmatrix(cc),j);
        %obj_img = imdilate(ismember(labelmatrix(cc),j),strel('disk',1));%'line',3,ss(j).Orientation));
        %obj_img = imerode(obj_img,strel('line',2,90+ss(j).Orientation));
        %obj_img = imdilate(obj_img,strel('line',ceil(0.05*(ss(j).MajorAxisLength)),ss(j).Orientation));
        rotated_object  = imrotate(obj_img,-angle);
        rotated_object = imdilate(rotated_object,strel('line',ceil(0.05*(ss(j).MajorAxisLength)),0));
        rotated_object = imerode(rotated_object,strel('line',ceil(0.05*(ss(j).MajorAxisLength)),0));
        rotated_object = imfill(rotated_object,'holes');
        filtered1 = imfilter(double(rotated_object),filter1,'conv');
        idx_down = find(filtered1 == -1);
        idx_up = find(filtered1 == 1);
        [idup_y,idup_x] = ind2sub(size(rotated_object),idx_up);
        [iddown_y,iddown_x] = ind2sub(size(rotated_object),idx_down);
        verticals = -idup_y+iddown_y;
        %newv =  medfilt1(verticals,3);
        newv = verticals;
        nwv = sort(newv);
         nwv1 = unique(nwv);
        %coef1(j) = mean(nwv(end-round(0.5*length(nwv)):end))/mean(newv);
         % coef1(j) = max(newv)/mean(newv);
         coef1(j) = nwv1(end-1)/mean(newv);
    end
    newiso = ismember(labelmatrix(cc),find(coef1<1.5));

    overlaps{i} = logical(overlaps{i}) + logical(isolated{i}) - newiso;
    isolatednew{i} = newiso;
    parallels{i} = logical(bw3-newiso);  
    
    bw_all = zeros(size(newiso));
  k = k+1;
    i1 = find(isolatednew{i} == 1);
    i1i(k) = length(i1);
    i2 = find(overlaps{i} == 1);
    i2o(k) = length(i2);
    i3 = find(border{i} == 1);
    i3b(k) = length(i3);
    %randomly assigned numbers to seperate them easily and make them colorful
     bw_all(i1) = uint8(250);
     bw_all(i2) = uint8(150);
     bw_all(i3) = uint8(75);
    bwrgb = uint8(bw_all);
    ll = multithresh(bwrgb,3); %segmentation is just for visual purposes
    bwrgb1 = imquantize(bwrgb,ll);
    bwrgb = label2rgb(bwrgb1);

    grouped{i} = bwrgb;
    
figure,subplot(121),imshow(image),subplot(122),imshow(bwrgb)

ans1 = questdlg('Apply the same grouping criteria for the rest of the images','Grouping check:','Yes','No','Yes');

k = 0;
if length(ans1) < 3
     return
else
    for i = 1:nfiles
        bw2 = pre_bw{i}; 
    s = regionprops(bw2, 'all');
    c = bwconncomp(bw2,8);
    for oi = 1:length(s)
        convimmini = s(oi).ConvexImage;
        convimmini1 = imrotate(convimmini,90-s(oi).Orientation);
        smini = regionprops(convimmini1,'BoundingBox');
        bboxvalues = smini.BoundingBox;
        mna(oi) = bboxvalues(3);
        mja(oi) = bboxvalues(4);
    end
    ar = mja./mna;
    bw3 = ismember(labelmatrix(c), find(ar>ar_min & mna < (max_minor/pixel_length) & mja > (min_major/pixel_length) & mja < (max_major/pixel_length)));
    bwoverlap1 = logical(bw2-bw3);
    bw3 = bwmorph(bw3,'thin',1);
    bwoverlap1 = bwareaopen(bwoverlap1,floor(numel(bwi)*0.0004));
    isolated{i} = bw3;
    overlaps{i} = (bwoverlap1);
    clusters = bwoverlap1;
    ss = regionprops(bw3,'all');
    cc = bwconncomp(bw3,8);
    image = ims{i};
    coef1 = [];
    coefnew = [];
    coefnewnew = [];
    coefnew_a = [];
    coefnew_b = [];
        
    for j = 1:length(ss)
        angle = ss(j).Orientation;
        obj_img = ismember(labelmatrix(cc),j);
        rotated_object  = imrotate(obj_img,-angle);
        rotated_object = imdilate(rotated_object,strel('line',ceil(0.05*(ss(j).MajorAxisLength)),0));
        rotated_object = imerode(rotated_object,strel('line',ceil(0.05*(ss(j).MajorAxisLength)),0));
        rotated_object = imfill(rotated_object,'holes');
        filtered1 = imfilter(double(rotated_object),filter1,'conv');
        idx_down = find(filtered1 == -1);
        idx_up = find(filtered1 == 1);
        [idup_y,idup_x] = ind2sub(size(rotated_object),idx_up);
        [iddown_y,iddown_x] = ind2sub(size(rotated_object),idx_down);
        verticals = -idup_y+iddown_y;
        newv =  verticals;
        nwv = sort(newv);
     nwv1 = unique(nwv);
        %coef1(j) = mean(nwv(end-round(0.5*length(nwv)):end))/mean(newv);
         % coef1(j) = max(newv)/mean(newv);
         coef1(j) = nwv1(end-1)/mean(newv);
    end
    newiso = ismember(labelmatrix(cc),find(coef1<1.5));

    overlaps{i} = logical(overlaps{i}) + logical(isolated{i}) - newiso;
    isolatednew{i} = newiso;
    parallels{i} = logical(bw3-newiso);  
    
  
    bw_all = ones(size(newiso));
    k = k+1;
    i1 = find(isolatednew{i} == 1);
    i1i(k) = length(i1);
    i2 = find(overlaps{i} == 1);
    i2o(k) = length(i2);
    i3 = find(border{i} == 1);
    i3b(k) = length(i3);
    %randomly assigned numbers to seperate them easily and make them colorful
     bw_all(i1) = uint8(250);
     bw_all(i2) = uint8(150);
     bw_all(i3) = uint8(75);
    bwrgb = uint8(bw_all);
    ll = multithresh(bwrgb,3); %segmentation is just for visual purposes
    bwrgb1 = imquantize(bwrgb,ll);
    bwrgb = label2rgb(bwrgb1);
  
    grouped{i} = bwrgb;
    
    f = figure('visible','off');
        subplot(121),imshow(image,[]),subplot(122),imshow(bwrgb)
        pname = path(i);
      name_ex = ['results\image_analysis_figures\', 'img',num2str(i),'grouped_duo''.png'];
     saveas(f,name_ex)    
    
    end
end

%%
i = 1;
k=0;
    bw2 = isolatednew{i}; 
    ss = regionprops(bw2, 'all');
    cc = bwconncomp(bw2,8);
    coef1 = [];
    lengthsup1 = [];
    widthsup25 = [];
    widthsup50 = [];
    l_max = [];
    w_max = [];

 f = figure('visible','on');
 imshow(ims{i},[]),hold on
    for j = 1:length(ss)
        angle = ss(j).Orientation;       
        obj_img = imdilate(ismember(labelmatrix(cc),j),strel('disk',1));%'line',3,ss(j).Orientation));   
        obj_img = imdilate(obj_img,strel('line',ceil(0.025*(ss(j).MajorAxisLength)),ss(j).Orientation));
        obj_img = imfill(obj_img,'holes');
        cell_coors = bwboundaries(obj_img);
        coors = cell_coors{1};
      cir2 =   coors(:,2) ;
      cir1 = coors(:,1);
     plot(cir2(1:3:end),cir1(1:3:end),'g-','LineWidth',1),hold on
      text(ss(j).Centroid(1),ss(j).Centroid(2),num2str(j),'Color','red','FontSize',15),hold on
    end
  allL = [];
  allW = [];
  headers = {'image_name','CNC_number','length (nm)','width (nm)','aspect ratio'};
  headers2 = {'image_name','isolated cnc area over total image area','agg. cnc area over total image area','border cnc area over total image area'};  
 cnc_numbers = [];
 image_names = {};
  for i = 1:nfiles
     f = figure('visible','off');
 imshow(ims{i},[]),hold on
 bw2 = isolatednew{i}; 
    ss = regionprops(bw2, 'all');
    cc = bwconncomp(bw2,8);
    lengthsup1 = [];
    widthsup25 = [];
    widthsup50 = [];
    l_max = [];
    w_max = [];
  for j = 1:length(ss)
        angle = ss(j).Orientation;       
        obj_img = imdilate(ismember(labelmatrix(cc),j),strel('disk',1));%'line',3,ss(j).Orientation));   
        obj_img = imdilate(obj_img,strel('line',ceil(0.025*(ss(j).MajorAxisLength)),ss(j).Orientation));
        obj_img = imfill(obj_img,'holes');
        cell_coors = bwboundaries(obj_img);
        coors = cell_coors{1};
      cir2 =   coors(:,2) ;
      cir1 = coors(:,1);
     plot(cir2(1:3:end),cir1(1:3:end),'g-','LineWidth',1),hold on
      text(ss(j).Centroid(1),ss(j).Centroid(2),num2str(j),'Color','red','FontSize',15),hold on
      rotated_object  = imrotate(obj_img,-angle);       
        rotated_object = imdilate(rotated_object,strel('line',ceil(0.05*(ss(j).MajorAxisLength)),0));
        rotated_object = imerode(rotated_object,strel('line',ceil(0.025*(ss(j).MajorAxisLength)),0));
        rotated_object = imfill(rotated_object,'holes');        
        filtered1 = imfilter(double(rotated_object),filter1,'conv');
        filtered2 = imfilter(double(rotated_object),filter2,'conv');
        idx_down = find(filtered1 == -1);
        idx_up = find(filtered1 == 1);
        [idup_y,idup_x] = ind2sub(size(rotated_object),idx_up);
        [iddown_y,iddown_x] = ind2sub(size(rotated_object),idx_down);
        idx_left = find(filtered2 == 1);
        idx_right = find(filtered2 == -1);
        [idleft_y,idleft_x] = ind2sub(size(rotated_object),idx_left);
        [idright_y,idright_x] = ind2sub(size(rotated_object),idx_right);
        [newlefty,Ilefty] = sort(idleft_y);
        newleftx = idleft_x(Ilefty);
        [newrighty,Irighty] = sort(idright_y);
        newrightx = idright_x(Irighty);
        verticals = -idup_y+iddown_y;
        horizontals = newrightx - newleftx;
        horizontals = medfilt1(horizontals,3);
        verticals = medfilt1(verticals,3);       
        lengthsup1(j) = pixel_length*mean(horizontals(horizontals>(max(horizontals)-0.1*range(horizontals))));
 widthsup50(j) = pixel_length*mean(verticals(verticals>(max(verticals)-0.50*range(verticals))));
 forname = struct2cell(dir(path{i}));
  image_names{(i-1)*length(ss) + j} = forname(1);
  end
    allL = [allL,lengthsup1];
    allW = [allW,widthsup50];
     name_ex = ['results\image_analysis_figures\', 'CNC_marked_img',num2str(i),'.png'];
     saveas(f,name_ex)      
     cnc_numbers = [cnc_numbers,1:j];
     %T2 = 
 end
 filename = 'results\measurements.xlsx';
aspect_ratio = [allL./allW]';
length_nm = allL';
width_nm = allW';
img_name = image_names';
cnc_no = cnc_numbers';
 T = table(img_name,cnc_no,length_nm,width_nm,aspect_ratio);
writetable(T,filename,'Sheet',1,'WriteRowNames',true)

%%
if length(allL) > 40
f = figure('visible','on');
subplot(131), histogram(allL,'Normalization','probability','BinWidth',10)
xlabel('length (nm)'),ylabel('probability (count/total N)'),ylim([0 0.16])
subplot(132), histogram(allW,'Normalization','probability','BinWidth',0.5)
xlabel('width (nm)'),ylabel('probability (count/total N)'),ylim([0 0.16])
subplot(133),histogram(allL./allW,20,'Normalization','probability')
xlabel('aspect ratio'),ylabel('probability (count/total N)')

for i = 1:length(allL)
    newL = [newL,allL(i)];
    meanL(i) = mean(newL);
    newW = [newW,allW(i)];
    meanW(i) = mean(newW);
end

allL2d = [allL,0,250];
    allW2d = [allW,0,15];
    figure
h1 = histogram2(allW2d,allL2d,'BinWidth',[0.5 10],'Normalization','probability','FaceColor','flat','ShowEmptyBins','on')
hv7 = h1.Values;
xlabel('width (nm)','FontSize',16)
ylabel('length (nm)','FontSize',16)
colormap jet
view(2)
axis([0 15 0 250])

figure,subplot(121),yyaxis left
plot(meanL,'*'),ylim([60 105]),ylabel('mean length (nm)'),hold on,plot(1:length(allL),ones(length(allL),1)*mean(allL),'k--','LineWidth',2.5), hold off
yyaxis right
plot(meanW,'*'),ylim([6 15]),hold on,xlim([0 length(allL)]),xlabel('number of measured CNCs')
ylabel('mean width (nm)'),plot(1:length(allW),ones(length(allW),1)*mean(allW),'k--','LineWidth',2.5)
hold off,grid on
subplot(122),plot(100*abs(meanL-mean(allL))./mean(allL),'*'),hold on,plot(100*abs(meanW-mean(allW))./mean(allW),'*')
xlabel('number of measured CNCs'),ylabel('percentage (%) change in mean'),ylim([0 15]),xlim([0 length(allW)])
legend('length','width'),grid on
end
%%







 
