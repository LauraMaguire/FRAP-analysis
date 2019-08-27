% Load images for a gel

%r = load('/Volumes/houghgrp/Processed Images/2019-2-6_3/results.mat');
r = load('/Volumes/houghgrp/Processed Images/2019-2-19_21/results.mat');% this is gel 34
%% Load the final results (does not have all images, that's why you need r too)
load('/Users/lauramaguire/Google Drive/Hough Lab/hydrogel-paper/data_final_190605.mat');
%% Set the index of the gel you care about
num=34;
%% Create a reconstruction at 5 min into the experiment (4 min after first bleach image)
% the 306 in the argument is time in minutes; corresponds to frame 19 for
% most experiments (closest possible to 300 s = 5 min), but some weird
% experiments might need manual adjusting here and of the frame number
% later
% Note that this does not take the 5 min experimental image as an input.
% These reconstructions only use the initial distribution then run time
% forward using the mode coefficients and diffusion constant from the fit.
rawRecon5min = calcTimeDist(im2double(data{num}.greenImage),data{num}.gelMask,...
    data{num}.cosArrayGrn,data{num}.sinArrayGrn,data{num}.numFitGrn2.D,306);
%% process the 5 min reconstruction
recon5min = rawRecon5min;
% I need to make the reservoir all NaNs so I can set a good scale using
% stretchlim, but I put the reservoir values back semicorrectly in the final
% images.  It's not the most efficient code.
recon5min(recon5min==0)=NaN;
recon5min = recon5min+data{num}.refGrn;
% make a good approximation of the background (reservoir) using the
% equilibrated intensity in the gel refGrn and the partition coefficient
% partC
background = data{num}.refGrn/data{num}.partC(1);
recon5min(isnan(recon5min))=background;
image5min = data{num}.gelMask.*im2double(r.GreenImages{19});%-data{34}.refGrn;
image5min(image5min==0) = NaN;

%% Create a reconstruction at about 1 hr into the experiment
rawRecon1hr = calcTimeDist(im2double(data{num}.greenImage),data{num}.gelMask,...
    data{num}.cosArrayGrn,data{num}.sinArrayGrn,data{num}.numFitGrn2.D,4026);
%% Process the 1 hr reconstruction
recon1hr = rawRecon1hr;
recon1hr(recon1hr==0)=NaN;
recon1hr = recon1hr+data{num}.refGrn;
background = data{num}.refGrn/data{num}.partC(1);
recon1hr(isnan(recon1hr))=background;
image1hr = data{num}.gelMask.*im2double(r.GreenImages{end});%-data{34}.refGrn;
image1hr(image1hr==0) = NaN;

%% Process the 1 min reconstruction (this is just the initialDistribution already
% saved in the data structure - first post-bleach image, taken about a
% minute after the end of the bleach)
rawRecon1min = data{num}.IDGrn;
recon1min = rawRecon1min;
recon1min(recon1min==0)=NaN;
recon1min = recon1min+data{num}.refGrn;
background = data{num}.refGrn/data{num}.partC(1);
recon1min(isnan(recon1min))=background;
image1min = data{num}.gelMask.*im2double(r.GreenImages{2});%-data{34}.refGrn;
image1min(image1min==0) = NaN;

%% Process the pre-bleach image
%imagePreBleach = data{num}.gelMask.*im2double(r.GreenImages{1});
imagePreBleach = im2double(r.GreenImages{1});
%imagePreBleach(imagePreBleach==0) = NaN;

%% Put everything into a figure
scale = 500/1.58; %1.58 um/px 4x Olympus
limIm = [min(min(image1min)); max(max(image1min))];
limRec = [min(min(recon1min)); max(max(recon1min))];
blank = zeros(1024,1344);

subplot(4,2,1)
%imshow(imadjust(imagePreBleach,limIm))
imshow(imadjust(cat(3,blank,imagePreBleach,blank),limIm))
line([100,100+scale],[100,100]) % show scale bar

subplot(4,2,3)
%imshow(imadjust(image1min,limIm))
%imshow(imadjust(im2double(r.GreenImages{2}),limIm))
imshow(imadjust(cat(3,blank,im2double(r.GreenImages{2}),blank),limIm))
%line([100,100+scale],[100,100]) % show scale bar

subplot(4,2,4)
%imshow(imadjust(recon1min,limRec))
imshow(imadjust(cat(3,blank,recon1min,blank),limIm))

subplot(4,2,5)
%imshow(imadjust(image5min,limIm))
%imshow(imadjust(im2double(r.GreenImages{19}),limIm))
imshow(imadjust(cat(3,blank,im2double(r.GreenImages{19}),blank),limIm))

subplot(4,2,6)
%imshow(imadjust(recon5min,limRec));
imshow(imadjust(cat(3,blank,recon5min,blank),limIm))

subplot(4,2,7)
%imshow(imadjust(image1hr,limIm))
%imshow(imadjust(im2double(r.GreenImages{end}),limIm))
imshow(imadjust(cat(3,blank,im2double(r.GreenImages{end}),blank),limIm))

subplot(4,2,8)
%imshow(imadjust(recon1hr,limRec))
imshow(imadjust(cat(3,blank,recon1hr,blank),limIm))

%% The following sections are for plotting the data and recovery fit together
% take the string that I fit to and make it a function
fg = inline(data{num}.KGrn2,'D','t','c1','c2');
fr = inline(data{num}.KRed2,'D','t','c1','c2');
%% create simulated data with the functions defined above
for t=1:length(data{num}.time)

dataSimTestG(t) = fg(data{num}.numFitGrn2.D,data{num}.time(t),...
data{num}.numFitGrn2.c1,data{num}.numFitGrn2.c2);
%dataSimTestG = dataSimTestG(:,1);
dataSimTestR(t) = fr(data{num}.numFitRed2.D,data{num}.time(t),...
data{num}.numFitRed2.c1,data{num}.numFitRed2.c2);
%dataSimTestR = dataSimTestR(:,1);

end
%% Plot things in nice colors for figures
plot(data{num}.time(2:end)/60,data{num}.norm(1,(2:end))/data{num}.norm(1,1),...
    'o','MarkerFaceColor',[41/255,116/255,81/255],...
    'MarkerEdgeColor',[41/255,116/255,81/255]);
hold on
plot(data{num}.time(2:end)/60,dataSimTestG(2:end)/data{num}.norm(1,1),...
    'color',[45/255,210/255,132/255]);
plot(data{num}.time(2:end)/60,data{num}.norm(2,(2:end))/data{num}.norm(2,1),...
    'o','MarkerEdgeColor',[235/255,45/255,45/255]);
plot(data{num}.time(2:end)/60,dataSimTestR(2:end)/data{num}.norm(2,1),...
    'color',[155/255,45/255,45/255]);
hold off
xlabel('Time (min)')
ylabel('Normalized intensity')
legend({'NTF2 data','NTF2 fit','mCherry data','mCherry fit'},'Location','southeast');
pbaspect([147.3 126 1])

