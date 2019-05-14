%% This and several following sections stitch together two segments of an experiment.
info = load('/Volumes/houghgrp/Processed Images/2019-1-18_1/2019-1-18_1--info.mat');
info2 = load('/Volumes/houghgrp/Processed Images/2019-1-18_2/2019-1-18_2--info.mat');
info = info.info;
info2 = info2.info;


% Adjust file paths based on PC/Mac difference.

if ispc % If this computer is a PC
    slash = '\'; % use backslashes along path
    expFolder = info.expFolderPC; % set base path to experiment
    baseSavePath = info.baseSavePathPC; % set base save path
else % If this computer is a Mac
    slash = '/'; % use forward slashes along path
    expFolder = info.expFolderMac; % set base path to experiment
    baseSavePath = info.baseSavePathMac; % set base save path
end

% Define useful variables from info structure.

greenConc = info.greenConc;
redConc = info.redConc;
greenName = info.greenName;
redName = info.redName;
conc = info.conc;
protein = info.protein;
geo = info.geo;
linker = info.linker;

% Import the video.

disp('Importing experiment...');
data = bfopen([expFolder slash info.expName]); % load experiment
info.frames = size(data{1,1},1)/info.nChannels; % verify number of frames
% save(path, 'info'); % overwrite old saved info structure with correct number of frames
% display(['Experiment has been imported. It has ', num2str(info.frames), ' frames.'])

disp('Importing experiment...');
data2 = bfopen([expFolder slash info2.expName]); % load experiment
info.frames2 = size(data2{1,1},1)/info2.nChannels; % verify number of frames

%%
% split images into green and red
ListOfImages1 = data{1,1};
GreenImages1 = ListOfImages1(1:2:length(ListOfImages1));
RedImages1 = ListOfImages1(2:2:length(ListOfImages1));

% split images into green and red
ListOfImages2 = data2{1,1};
GreenImages2 = ListOfImages2(1:2:length(ListOfImages2));
RedImages2 = ListOfImages2(2:2:length(ListOfImages2));

% % define the contrast scale based on the first frame and carry that scale
% % through the entire movie
grnScale = stretchlim(im2double(GreenImages1{1,1}));
redScale = stretchlim(im2double(RedImages1{1,1}));

RedImages = [RedImages1 RedImages2];
GreenImages = [GreenImages1 GreenImages2];

%%
sav = '/Volumes/houghgrp/Processed Images/2019-1-18_1';
imshow(imadjust(im2double(RedImages{1,1}),redScale));
[~, ~, total, ~, ~] = roipoly(imadjust(im2double(RedImages{1,1}),redScale));
imshow(total);
save([sav 'total'],'total');
close all

%%
%time = [5*(1:60) 300+60*(1:10)];
%time = 5*(1:180);
segment1 = 18+5*(1:120);
segment2 = 58+20*(1:20);
time = [segment1 (segment1(end)+segment2)];
%%
n = 140;

%%
totr = zeros(n,1);
totg = zeros(n,1);
bleach = zeros(n,1);
bleach2 = zeros(n,1);
for t=1:n
    totr(t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*total));
    totg(t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*total));
    areaRatio=(sum(sum(total))/sum(sum(bleachSpot)));
    bleach(t) = areaRatio*sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*bleachSpot))/totr(t);
    bleach2(t) = areaRatio*sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*bleachSpot))/totg(t);
    disp(t);
end

plot(time,totr/totr(1),'ro');
hold on
plot(time,totg/totg(1),'go');
hold off

figure
plot(time,bleach,'ro');
hold on
plot(time,bleach2,'go');
hold off

%%
[fitresult, gof] = FRAPfit(time, bleach);
a = fitresult.a;
c = fitresult.c;
tau = fitresult.tau;

%%

%%
bleach = zeros(n,1);
ref = zeros(n,1);
for t=1:n
    bleach(t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*bleachSpot))/sum(sum(bleachSpot));
    ref(t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*refSpot))/sum(sum(refSpot));
    disp(t);
end

bleach2 = zeros(n,1);
ref2 = zeros(n,1);
for t=1:n
    bleach2(t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*bleachSpot))/sum(sum(bleachSpot));
    ref2(t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*refSpot))/sum(sum(refSpot));
    disp(t);
end


%%
sav = '/Volumes/houghgrp/Processed Images/2019-1-3_12/';

%%
mov = VideoWriter([sav 'FRAP.avi']);
% Set the frame rate. Always will make a 10-s movie.
mov.FrameRate = info.frames/10;

% Open the video writer before beginning.
open(mov);

% Loop over all green images, display them, and add them to the movie.
for i=1:70
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
    disp(['Frame ' num2str(i) ' of ' num2str(info.frames)]);
end

%Close the final figure.
close all

%Close the video writer.
close(mov);


%%
%%
sav = '/Volumes/houghgrp/Processed Images/2019-1-3_11/';
imshow(imadjust(im2double(RedImages{1,1}),redScale));
[~, ~, total, ~, ~] = roipoly(imadjust(im2double(RedImages{1,1}),redScale));
imshow(total);
save([sav 'total'],'total');

%%
n = 60;
totr = zeros(n,1);
totg = zeros(n,1);
for t=1:n
    totr(t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*total));
    totg(t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*total));
    disp(t);
end

plot(totr/totr(1),'r');
hold on
plot(totg/totg(1),'g');
hold off
%%
imshow(imadjust(im2double(RedImages{1,1}),redScale));
[~, ~, bleachSpot, ~, ~] = roipoly(imadjust(im2double(RedImages{1,1}),redScale));

imshow(bleachSpot);
save([sav 'bleachSpot'],'bleachSpot');
close all

%%
imshow(imadjust(im2double(RedImages{1,1}),redScale));
[~, ~, refSpot, ~, ~] = roipoly(imadjust(im2double(RedImages{1,1}),redScale));
%%
imshow(refSpot);
save([sav 'refSpot'],'refSpot');

%%
bleach = zeros(60,1);
ref = zeros(60,1);
for t=1:60
    bleach(t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*bleachSpot))/sum(sum(bleachSpot));
    ref(t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*refSpot))/sum(sum(refSpot));
    disp(t);
end

%%
bleach2 = zeros(60,1);
ref2 = zeros(60,1);
for t=1:60
    bleach2(t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*bleachSpot))/sum(sum(bleachSpot));
    ref2(t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*refSpot))/sum(sum(refSpot));
    disp(t);
end

%%
%time = [5*(1:60) 303+60*(1:10)];
time = 5*(1:60);
norm = bleach2./test2;% bleach./ref;
plot(time,norm,'o');


%%
plot(bleach/bleach(end));
hold all

plot(ref)
plot(bleach2/bleach2(end))
plot(ref2)

%%
w = (1.58)*sqrt(sum(sum(bleachSpot))/pi)
t12 = -log(0.5)*529
D = 0.88*w^2/(4*t12)
%%
surf(imadjust(im2double(RedImages{1,1}),redScale));
%%
mesh(1:1344,1:1024,imadjust(im2double(RedImages{1,1}),redScale));

%%
x=1:474;
y=1:804;

%%
J = imcrop(im2double(RedImages{1,1}),redScale);

%%
test2 = sum(sum(bleachSpot))*bleach2+sum(sum(refSpot))*ref2;
%%
plot(test2);
hold all
plot(test);

%%
sum(sum(bleachSpot))/(sum(sum(bleachSpot))+sum(sum(refSpot)));