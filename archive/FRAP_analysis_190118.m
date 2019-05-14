%% This and several following sections stitch together two segments of an experiment.
info1 = load('/Volumes/houghgrp/Processed Images/2019-2-28_15/2019-2-28_15--info.mat');
info2 = load('/Volumes/houghgrp/Processed Images/2019-2-28_16/2019-2-28_16--info.mat');
sav = '/Volumes/houghgrp/Processed Images/2019-2-28_15/';

prebleach = '/Volumes/houghgrp/Microscopy/190228_FRAP/190228_FRAP_26.vsi';

offset1 = 35; % time between end of bleaching and start of recovery sequence, in seconds
nFrames1 = 30; % number of frames in segment 1
deltaT1 = 15; % time per frame in segment 1, in seconds

offset2 = 60; % time between end of seg. 1 and start of seg. 2, in seconds
nFrames2 = 30; % number of frames in segment 2
deltaT2 = 120; %time per frame in segment 2, in seconds

%%
info1 = info1.info;
info2 = info2.info;

segmentpre = [1];

segment1 = offset1+deltaT1*(1:nFrames1);
segment2 = offset2+deltaT2*(0:nFrames2-1);

time = [segmentpre (segmentpre+segment1) (segmentpre+segment1(end)+segment2)];

%time = [segment1 (segment1(end)+segment2)];

%time = [5*(1:60) 300+60*(1:10)];
%time = 10+5*(1:180);
% segment1 = 18+10*(1:30);
% segment2 = 58+60*(0:19);
% time = [segment1 (segment1(end)+segment2)];

n = length(time);

clear segment1 segment2

% Adjust file paths based on PC/Mac difference.

if ispc % If this computer is a PC
    slash = '\'; % use backslashes along path
    expFolder = info1.expFolderPC; % set base path to experiment
    baseSavePath = info1.baseSavePathPC; % set base save path
else % If this computer is a Mac
    slash = '/'; % use forward slashes along path
    expFolder = info1.expFolderMac; % set base path to experiment
    baseSavePath = info1.baseSavePathMac; % set base save path
end

% Define useful variables from info structure.

greenConc = info1.greenConc;
redConc = info1.redConc;
greenName = info1.greenName;
redName = info1.redName;
conc = info1.conc;
protein = info1.protein;
geo = info1.geo;
linker = info1.linker;

%% Import the video.

disp('Importing experiment...');
data = bfopen([expFolder slash info1.expName]); % load experiment
info1.frames = size(data{1,1},1)/info1.nChannels; % verify number of frames
% save(path, 'info'); % overwrite old saved info structure with correct number of frames
% display(['Experiment has been imported. It has ', num2str(info.frames), ' frames.'])

disp('Importing experiment...');
data2 = bfopen([expFolder slash info2.expName]); % load experiment
info1.frames2 = size(data2{1,1},1)/info2.nChannels; % verify number of frames

% disp('Importing experiment...');
datapre = bfopen(prebleach); % load experiment

%%
% split images into green and red
ListOfImagesPre = datapre{1,1};
GreenImagespre = ListOfImagesPre(1:2:length(ListOfImagesPre));
RedImagespre = ListOfImagesPre(2:2:length(ListOfImagesPre));

ListOfImages1 = data{1,1};
GreenImages1 = ListOfImages1(1:2:length(ListOfImages1));
RedImages1 = ListOfImages1(2:2:length(ListOfImages1));

% split images into green and red
ListOfImages2 = data2{1,1};
GreenImages2 = ListOfImages2(1:2:length(ListOfImages2));
RedImages2 = ListOfImages2(2:2:length(ListOfImages2));

% % define the contrast scale based on the first frame and carry that scale
% % through the entire movie
grnScale = stretchlim(im2double(GreenImages1{1,2}));
redScale = stretchlim(im2double(RedImages1{1,2}));

RedImages = [RedImagespre RedImages1 RedImages2];
GreenImages = [GreenImagespre GreenImages1 GreenImages2];
% RedImages = [RedImages1 RedImages2];
% GreenImages = [GreenImages1 GreenImages2];
% RedImages = [RedImages1];
% GreenImages = [GreenImages1];
%%
mov = VideoWriter([sav 'FRAP.avi']);
% Set the frame rate. Always will make a 10-s movie.
mov.FrameRate = info1.frames/10;

% Open the video writer before beginning.
open(mov);

% Loop over all green images, display them, and add them to the movie.
for i=1:n
    % Convert to from 16-bit image data to double.
    imageg = imadjust(im2double(GreenImages{1,i}),grnScale);
    imager = imadjust(im2double(RedImages{1,i}),redScale);
    imageb = zeros(size(imager));
    image = cat(3,imager,imageg,imageb);
    % Display the image as a figure.
    imshow(image, 'InitialMagnification',50);
    % Make the figure a movie frame. ('gcf' means 'get current figure')
    F = getframe(gcf);
    % Add the frame to the movie.
    writeVideo(mov,F);
    disp(['Frame ' num2str(i)]);
end

%Close the final figure.
close all

%Close the video writer.
close(mov);

%% mask the whole gel
imshow(imadjust(im2double(RedImages{1,2}),redScale));
[~, ~, total, ~, ~] = roipoly(imadjust(im2double(RedImages{1,2}),redScale));
imshow(total);
save([sav 'total'],'total');
close all

%% mask just the bleach spot
imshow(imadjust(im2double(RedImages{1,2}),redScale));
[~, ~, bleachSpot, ~, ~] = roipoly(imadjust(im2double(RedImages{1,2}),redScale));
imshow(bleachSpot);
save([sav 'bleachSpot'],'bleachSpot');
close all

%%
finalGreen = im2double(GreenImages{1,end});
finalGreen = imadjust(finalGreen,grnScale);

finalRed = im2double(RedImages{1,2});
finalRed = imadjust(finalRed,redScale);

composite = imfuse(finalGreen, finalRed, 'falsecolor');
imshow(composite)
%%
[~, ~, bleachSpot, ~, ~] = roipoly(composite);
imshow(bleachSpot);
save([sav 'bleachSpot'],'bleachSpot');
close all
%% mask the gel minus the bleach spot

imshow(imadjust(im2double(RedImages{1,2}),redScale));
[~, ~, refSpot, ~, ~] = roipoly(imadjust(im2double(RedImages{1,2}),redScale));
imshow(refSpot);
save([sav 'refSpot'],'refSpot');
close all


