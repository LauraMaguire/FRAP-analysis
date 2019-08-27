%% dataset used in thesis (1:43 actually used, above 43 not used)

folders = {'19-1-3_1','19-1-3_3','19-1-3_6','19-2-6_1','19-2-6_3',...
    '19-2-6_5','19-2-6_7','19-2-6_9','19-2-6_11','19-2-6_13',...
    '19-2-6_15','19-2-6_17','19-2-6_19','19-2-6_21','19-2-6_23',...
    '19-2-12_3','19-2-12_5','19-2-12_7','19-2-12_9','19-2-12_11',...
    '19-2-12_13','19-2-12_15','19-2-12_17','19-2-15_1','19-2-15_7',...
    '19-2-15_13','19-2-15_15','19-2-15_19','19-2-15_21','19-2-19_1',...
    '19-2-19_7','19-2-19_15','19-2-19_19','19-2-19_21','19-2-19_23',...
    '19-2-21_7','19-2-21_11','19-2-21_13','19-2-21_15','19-2-21_17',...
    '19-2-21_19','19-2-21_21','19-2-21_23','19-2-28_1','19-2-28_3'...
    '19-2-28_5','19-2-28_7','19-2-28_9','19-2-28_11','19-2-28_13','19-2-28_15'}.';
data = cell(length(folders),1);

%% test dataset of high-concentration gels (2-28 and 3-1)
% foldersTest = {'19-2-28_1','19-2-28_3','19-2-28_5','19-2-28_7','19-2-28_9',...
%     '19-2-28_11','19-2-28_13','19-2-28_15','19-3-1_1','19-3-1_3','19-3-1_5'}.';
% data = cell(length(foldersTest),1);

%% dataset in thesis (1:43)
type = {'ctrl','cct1','ctrl','ctrl','cct1','cct2','ctrl','cct1','cct2',...
    'ctrl','cct1','cct2','ctrl','cct1','cct2','cct1','cct2','ctrl','cct1',...
    'ctrl','cct1','cct2','cct2','ctrl','ctrl','ctrl','cct1','ctrl','cct1',...
    'ctrl','ctrl','cct1','ctrl','cct1','cct2','ctrl','cct2','ctrl','cct1',...
    'cct2','ctrl','cct1','cct2','ctrl','cct1','cct2','ctrl','cct2','ctrl','cct1','cct2'};
conc = [0 10 0 0 10 10 0 10 10 0 10 10 0 10 10 10 10 0 10 0 10 10 10 0 0 0 ...
    10 0 10 0 0 10 0 10 10 0 5 0 10 5 0 10 5 0 30 40 0 40 0 30 40];

%% test dataset
% type = {'ctrl','','cct1','cct2','ctrl','cct2','ctrl','cct1','cct2',...
%     'ctrl','cct1','cct2'};
% conc = [0 30 40 0 40 0 30 40 0 30 40];
%% run this to make the data structure needed for Fourier series fit
% The gel masks are, in order: entire gel, bleach spot, reservoir,
% equilibrated areas of gel only.
for n=1:length(data)
    % load 'results' structure that was the result of runnin
    % FRAP_processing_190513; takes a few minutes
    r = load(['/Volumes/houghgrp/Processed Images/20' foldersTest{n} '/results.mat']);
    
    % set scales for viewing images, NOT for calculating recovery curves
    data{n}.grnScale = stretchlim(im2double(r.GreenImages{2}));
    data{n}.redScale = stretchlim(im2double(r.RedImages{2}));
    
    % take first post-bleach image and save in data structure
    data{n}.greenImage = r.GreenImages{2};
    greenImage = im2double(r.GreenImages{2});
    data{n}.redImage = r.RedImages{2};
    redImage = im2double(r.RedImages{2});
    
    % make adjusted images for viewing ONLY
    greenImageAdj = imadjust(greenImage,data{n}.grnScale);
    redImageAdj = imadjust(redImage,data{n}.redScale);
    blueImage = zeros(size(redImage));
    composite = cat(3,redImageAdj,greenImageAdj,blueImage);
    
    % Mask the entire gel
    %imshow(imadjust(im2double(r.RedImages{1,2}),data{n}.redScale));
    [~, ~, gelMask, ~, ~] = roipoly(composite);
    imshow(gelMask);
    data{n}.gelMask = gelMask;
    close all

    % Mask the bleach spot
    %imshow(composite)
    [~, ~, bleachSpot, ~, ~] = roipoly(composite);
    imshow(bleachSpot);
    data{n}.bleachSpot = bleachSpot;
    close all
    
    % Mask the reservoir
    %imshow(imadjust(im2double(r.RedImages{1,1}),data{n}.redScale));
    [~, ~, resRef, ~, ~] = roipoly(composite);
    imshow(resRef);
    data{n}.resRef = resRef;
    close all
    
    % Mask the parts of the gel which are well-equilibrated
    [~, ~, equilMask, ~, ~] = roipoly(composite);
    imshow(equilMask);
    data{n}.equilMask = equilMask;
    close all

    % Calculate the recovery curves
    bleachRecovery = zeros(2,n);
    wholeGel = zeros(2,n);

    for t=1:length(r.time)
        bleachRecovery(1,t) = sum(sum(im2double(r.GreenImages{1,t}).*data{n}.bleachSpot));
        bleachRecovery(2,t) = sum(sum(im2double(r.RedImages{1,t}).*data{n}.bleachSpot));

        wholeGel(1,t) = sum(sum(im2double(r.GreenImages{1,t}).*data{n}.gelMask));
        wholeGel(2,t) = sum(sum(im2double(r.RedImages{1,t}).*data{n}.gelMask));
        
        data{n}.bleachRecovery = bleachRecovery;
        data{n}.wholeGel = wholeGel;
    end
    
    % Ratios for normalizing recovery curves
    bleachWholeRatio = sum(sum(bleachSpot))/sum(sum(gelMask));
    data{n}.norm = (bleachRecovery./wholeGel)/bleachWholeRatio;
    
    % partition coefficient calculations, background values hard-coded
    % because they were very stable over a wide range of exposure time and
    % gain conditions (all at 100% laser power)
    totPart = sum(sum(im2double(data{n}.greenImage-214).*data{n}.equilMask));
    totRes = sum(sum(im2double(data{n}.greenImage-214).*data{n}.resRef));
    areaPart = sum(sum(data{n}.equilMask));
    areaRes = sum(sum(data{n}.resRef));
    data{n}.partCeq(1) = (totPart/areaPart)/(totRes/areaRes);
    
    totPart = sum(sum(im2double(data{n}.redImage-219).*data{n}.equilMask));
    totRes = sum(sum(im2double(data{n}.redImage-219).*data{n}.resRef));
    areaPart = sum(sum(data{n}.equilMask));
    areaRes = sum(sum(data{n}.resRef));
    data{n}.partCeq(2) = (totPart/areaPart)/(totRes/areaRes);
    
    % bound probability calculations
    data{n}.bProbeq = 1-(data{n}.partCeq(2)/data{n}.partCeq(1));
    data{n}.time = r.time;
    
    % reference values = average value of equilibrated portions of
    % hydrogel, background values hard-coded and subtracted
    data{n}.refGrn = sum(sum(im2double(data{n}.greenImage-214).*data{n}.equilMask))...
    /sum(sum(equilMask));
    data{n}.refRed = sum(sum(im2double(data{n}.redImage-219).*data{n}.equilMask))...
    /sum(sum(equilMask));

    data{n}.type = type{n};
    data{n}.conc = conc(n);
    
    disp(['Finished ' num2str(n) ' of ' num2str(length(data)) '.']);
end

clear bleachSpot composite finalGreen finalRed n t gelMask bleachRecovery...
    wholeGel bleachWholeRatio resRef areaPart totRes areaRes totPart...
    greenImageAdj blueImage equilMask greenImage r redImage redImageAdj
