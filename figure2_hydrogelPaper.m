%% profile info and plots structure
load('/Volumes/houghgrp/Processed Images/2019-2-20_2/2019-2-20_2--info.mat');
load('/Volumes/houghgrp/Processed Images/2019-2-20_2/2019-2-20_2--plots.mat');
%% FRAP results structure
r = load('/Volumes/houghgrp/Processed Images/2019-2-21_5/results.mat');
%% Raw images
images = bfopen('/Volumes/houghgrp/Microscopy/190220/190220_profiles_02.vsi');
%%
images = images{1,1};
imagesGrn = images(1:2:length(images));
imagesRed = images(2:2:length(images));
clear images
%%

d.imageInitGrn = im2double(imagesGrn{1});
d.imageInitRed = im2double(imagesRed{1});
d.imageFinlGrn = im2double(r.GreenImages{1});
d.imageFinlRed = im2double(r.RedImages{1});

%% Define the gel mask for the profile experiment

grnScale = stretchlim(d.imageInitGrn);
redScale = stretchlim(d.imageInitRed);
    
greenImageAdj = imadjust(d.imageInitGrn,grnScale);
redImageAdj = imadjust(d.imageInitRed,redScale);
blueImage = zeros(size(d.imageInitRed));
composite = cat(3,redImageAdj,greenImageAdj,blueImage);

[~, ~, d.gelMask, ~, ~] = roipoly(composite);
imshow(d.gelMask);
close all
clear grnScale redScale greenImageAdj redImageAdj blueImage composite

%% Find the equilibrated concentration using the prebleach FRAP images
d.refGrn = mean(mean(d.gelMask.*d.imageFinlGrn));
d.refRed = mean(mean(d.gelMask.*d.imageFinlRed));

%% Fill in some more stuff
d.timeAx = plots.timeAx;

d.normAccGrn = plots.accumulation(1,:)./plots.reservoir(1,:);
d.normAccRed = plots.accumulation(2,:)./plots.reservoir(2,:);

d.accumulation = plots.accumulation;
d.reservoir = plots.reservoir;

[d.x,d.y] = findCenter(d.gelMask);

%%
n = 20;
[d.cosArrayGrn, d.sinArrayGrn, d.rmax] = calculateCoeffs(d.imageInitGrn-d.refGrn, d.gelMask, n, n);
[d.cosArrayRed, d.sinArrayRed, d.rmax] = calculateCoeffs(d.imageInitRed-d.refRed, d.gelMask, n, n);
d.IDGrn = calcInitDist(d.imageInitGrn, d.gelMask, d.cosArrayGrn, d.sinArrayGrn);
d.IDRed = calcInitDist(d.imageInitRed, d.gelMask, d.cosArrayRed, d.sinArrayRed);
%%
d.fitStringGrn = generateFitString(d.imageInitGrn, d.gelMask, ...
    d.cosArrayGrn, d.sinArrayGrn, d.rmax, d.x, d.y);
d.normFitStringGrn = ['c1*(' fitStringGrn ')+c2'];

[d.fitresultGrn, d.gofGrn] = numericalBesselFitAcc(d, d.normFitStringGrn,1);
%%
d.fitStringRed = generateFitString(d.imageInitRed, d.gelMask, ...
    d.cosArrayRed, d.sinArrayRed, d.rmax, d.x, d.y);
%%
d.normFitStringRed = ['c1*(' d.fitStringRed ')+c2'];

[d.fitresultRed, d.gofRed] = numericalBesselFitAcc(d, d.normFitStringRed,2);
%% next do the FRAP analysis - mostly cut and paste from other code
% set scales for viewing images, NOT for calculating recovery curves
f.grnScale = stretchlim(im2double(r.GreenImages{2}));
f.redScale = stretchlim(im2double(r.RedImages{2}));
    
% take first post-bleach image and save in data structure
f.greenImage = r.GreenImages{2};
greenImage = im2double(r.GreenImages{2});
f.redImage = r.RedImages{2};
redImage = im2double(r.RedImages{2});
    
% make adjusted images for viewing ONLY
greenImageAdj = imadjust(greenImage,f.grnScale);
redImageAdj = imadjust(redImage,f.redScale);
blueImage = zeros(size(redImage));
composite = cat(3,redImageAdj,greenImageAdj,blueImage);
    
% Mask the entire gel
[~, ~, gelMask, ~, ~] = roipoly(composite);
imshow(gelMask);
f.gelMask = gelMask;
close all

% Mask the bleach spot
[~, ~, bleachSpot, ~, ~] = roipoly(composite);
imshow(bleachSpot);
f.bleachSpot = bleachSpot;
close all
    
% Mask the reservoir
[~, ~, resRef, ~, ~] = roipoly(composite);
imshow(resRef);
f.resRef = resRef;
close all
    
% Mask the parts of the gel which are well-equilibrated
[~, ~, equilMask, ~, ~] = roipoly(composite);
imshow(equilMask);
f.equilMask = equilMask;
close all

% Calculate the recovery curves
bleachRecovery = zeros(2,n);
wholeGel = zeros(2,n);

for t=1:length(r.time)
    bleachRecovery(1,t) = sum(sum(im2double(r.GreenImages{1,t}).*f.bleachSpot));
    bleachRecovery(2,t) = sum(sum(im2double(r.RedImages{1,t}).*f.bleachSpot));

    wholeGel(1,t) = sum(sum(im2double(r.GreenImages{1,t}).*f.gelMask));
    wholeGel(2,t) = sum(sum(im2double(r.RedImages{1,t}).*f.gelMask));
        
    f.bleachRecovery = bleachRecovery;
    f.wholeGel = wholeGel;
end
    
% Ratios for normalizing recovery curves
bleachWholeRatio = sum(sum(bleachSpot))/sum(sum(gelMask));
f.norm = (bleachRecovery./wholeGel)/bleachWholeRatio;
    
% partition coefficient calculations, background values hard-coded
% because they were very stable over a wide range of exposure time and
% gain conditions (all at 100% laser power)
totPart = sum(sum(im2double(f.greenImage-214).*f.equilMask));
totRes = sum(sum(im2double(f.greenImage-214).*f.resRef));
areaPart = sum(sum(f.equilMask));
areaRes = sum(sum(f.resRef));
f.partCeq(1) = (totPart/areaPart)/(totRes/areaRes);
    
totPart = sum(sum(im2double(f.redImage-219).*f.equilMask));
totRes = sum(sum(im2double(f.redImage-219).*f.resRef));
areaPart = sum(sum(f.equilMask));
areaRes = sum(sum(f.resRef));
f.partCeq(2) = (totPart/areaPart)/(totRes/areaRes);
    
% bound probability calculations
f.bProbeq = 1-(f.partCeq(2)/f.partCeq(1));
f.time = r.time;
    
% reference values = average value of equilibrated portions of
% hydrogel, background values hard-coded and subtracted
f.refGrn = sum(sum(im2double(f.greenImage-214).*f.equilMask))...
    /sum(sum(equilMask));
f.refRed = sum(sum(im2double(f.redImage-219).*f.equilMask))...
    /sum(sum(equilMask));
    
clear bleachSpot composite finalGreen finalRed t gelMask bleachRecovery...
    wholeGel bleachWholeRatio resRef areaPart totRes areaRes totPart...
    greenImageAdj blueImage equilMask greenImage redImage redImageAdj

%%
[f.x,f.y] = findCenter(f.gelMask);
    
wholeMask = f.gelMask;
bleachMask = f.bleachSpot;
refMask = wholeMask-bleachMask;
    
image = im2double(f.greenImage);
image2 = image-f.refGrn;
    
[cosArrayGrn, sinArrayGrn, ~] = calculateCoeffs(image2, wholeMask, n,n);
f.cosArrayGrn = cosArrayGrn;
f.sinArrayGrn = sinArrayGrn;
disp('Finished green coeffs');
     
image = im2double(f.redImage);
image2 = image-f.refRed;
    
[cosArrayRed, sinArrayRed, rmax] = calculateCoeffs(image2, wholeMask, n, n);
f.cosArrayRed = cosArrayRed;
f.sinArrayRed = sinArrayRed;
f.rmax = rmax;
disp('Finished red coeffs');
    
initialDistribution = calcInitDist(image, wholeMask, cosArrayGrn, sinArrayGrn);
f.IDGrn = initialDistribution;
disp('Finished green init. dist.');
    
initialDistribution = calcInitDist(image, wholeMask, cosArrayRed, sinArrayRed);
f.IDRed = initialDistribution;
disp('Finished red init. dist.');

cosArrayRed=f.cosArrayRed;
sinArrayRed=f.sinArrayRed;
cosArrayGrn=f.cosArrayGrn;
sinArrayGrn=f.sinArrayGrn;
    
f.fitStringSpotGrn = generateFitString(image, f.bleachSpot,...
        cosArrayGrn, sinArrayGrn, f.rmax, f.x, f.y);
    
f.fitStringSpotRed = generateFitString(image, f.bleachSpot,...
        cosArrayRed, sinArrayRed, f.rmax, f.x, f.y);
    
f.fitStringGelGrn = generateFitString(image, f.gelMask,...
        cosArrayGrn, sinArrayGrn, f.rmax, f.x, f.y);
    
f.fitStringGelRed = generateFitString(image, f.gelMask,...
        cosArrayRed, sinArrayRed, f.rmax, f.x, f.y);

% insert extra fit parameters to fit to bleach depth and asymptotic
% recovered value
f.fGrn = ['(' f.fitStringSpotGrn '/'...
        num2str(sum(sum(f.bleachSpot))) '+' num2str(f.refGrn)...
        ')/(' f.fitStringGelGrn '/' num2str(sum(sum(f.gelMask)))...
        '+' num2str(f.refGrn) ')'];
    
f.KGrn = ['c1*(' f.fGrn ')+c2'];
    
f.fRed = ['(' f.fitStringSpotRed '/'...
        num2str(sum(sum(f.bleachSpot))) '+' num2str(f.refRed)...
        ')/(' f.fitStringGelRed '/' num2str(sum(sum(f.gelMask)))...
        '+' num2str(f.refRed) ')'];
    
f.KRed = ['c1*(' f.fRed ')+c2'];

[~,f.numFitGrn, ~] = numericalBesselFit(f, f.KGrn,1);

[~,f.numFitRed, ~] = numericalBesselFit(f, f.KRed,2);


