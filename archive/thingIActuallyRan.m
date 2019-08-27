folders = {'/Volumes/houghgrp/Processed Images/2018-2-2_2/2018-2-2_2--info.mat',...
    '/Volumes/houghgrp/Processed Images/2018-5-16_3/2018-5-16_3--info.mat',...
    '/Volumes/houghgrp/Processed Images/2018-7-3_1/2018-7-3_1--info.mat',...
    '/Volumes/houghgrp/Processed Images/2018-7-3_2/2018-7-3_2--info.mat',...
    };

foldersP = {'/Volumes/houghgrp/Processed Images/2018-2-2_2/2018-2-2_2--plots.mat',...
    '/Volumes/houghgrp/Processed Images/2018-5-16_3/2018-5-16_3--plots.mat',...
    '/Volumes/houghgrp/Processed Images/2018-7-3_1/2018-7-3_1--plots.mat',...
    '/Volumes/houghgrp/Processed Images/2018-7-3_2/2018-7-3_2--plots.mat',...
    };

for n=1:length(folders)
    data{n}.info = load(folders{n},'info');
    data{n}.plots = load(foldersP{n},'plots');
end

%%
[recoveryCurveGel]=simulateDataLargeMask(data{39}.greenImage,data{39}.gelMask,...
    data{39}.cosArrayGrn,data{39}.sinArrayGrn,data{39}.rmax,...
    data{39}.time,data{39}.numFitGrn.D,data{39}.x,data{39}.y);
%%
plot(data{39}.time,recoveryCurve/sum(sum(data{39}.bleachSpot)))
figure
plot(data{39}.time,recoveryCurveGel/sum(sum(data{39}.gelMask)))
figure
plot(data{39}.time(2:end), data{39}.wholeGel(1,2:end),'o');
figure
imagesc(data{39}.IDGrn)
figure
imagesc(data{39}.greenImage)
%%
f = (recoveryCurve/sum(sum(data{39}.bleachSpot))+data{39}.refGrn)...
    ./(recoveryCurveGel/sum(sum(data{39}.gelMask))+data{39}.refGrn);
plot(data{39}.time,f)

%%
for n=1:43
data{n}.fGrn2 = ['(' data{n}.fitStringSpotGrn '/'...
    num2str(sum(sum(data{n}.bleachSpot))) '+' num2str(data{n}.refGrn)...
    ')/(' data{n}.fitStringGelGrn '/' num2str(sum(sum(data{n}.gelMask)))...
    '+' num2str(data{n}.refGrn) ')'];
data{n}.KGrn2 = ['c1*(' data{n}.fGrn2 ')+c2'];
data{n}.fRed2 = ['(' data{n}.fitStringSpotRed '/'...
    num2str(sum(sum(data{n}.bleachSpot))) '+' num2str(data{n}.refRed)...
    ')/(' data{n}.fitStringGelRed '/' num2str(sum(sum(data{n}.gelMask)))...
    '+' num2str(data{n}.refRed) ')'];
data{n}.KRed2 = ['c1*(' data{n}.fRed2 ')+c2'];
disp(n);
end
%%
good = [1 2 4:9 11:14 16:28 30:35 37 39:41 43];
%%
%for n=[5 7]
%for n=1:43
for n=[1 2 4:9 11:14 16:28 30:35 37 39:41 43]
[~,data{n}.numFitGrn2, ~] = numericalBesselFit(data{n}, data{n}.KGrn2,1);
end
%%
%for n=[5 7]
for n=[1 2 4:9 11:14 16:28 30:35 37 39:41 43]
[~,data{n}.numFitRed2, ~] = numericalBesselFit(data{n}, data{n}.KRed2,2);
end

%%
for n=1:43
    D2(n,1) = data{n}.numFitGrn2.D;
    D2(n,2) = data{n}.numFitRed2.D;
    c1(n,1) = data{n}.numFitGrn2.c1;
    c2(n,1) = data{n}.numFitGrn2.c2;
    c1(n,2) = data{n}.numFitRed2.c1;
    c2(n,2) = data{n}.numFitRed2.c2;
    Dold(n,1) = data{n}.D(1);
end
%%
plot(D2(indexCct1),'o')
hold on
plot(D2(indexCct1),'o')
plot(D2(indexCct2),'o')
hold off
%%
for i=1:43
    figure
    imagesc(data{i}.IDGrn)
end
%%
n=16;
for v=1:4
    vv=[2,4,10,20];
    [data{n}.x,data{n}.y] = findCenter(data{n}.gelMask);
    
    wholeMask = data{n}.gelMask;
    bleachMask = data{n}.bleachSpot;
    refMask = wholeMask-bleachMask;
    
    image = im2double(data{n}.greenImage);
    data{n}.refGrn = sum(sum(image.*refMask))/sum(sum(refMask));
    image2 = image-data{n}.refGrn;
    
    image2 = image-ref;
    [r{v}.cosArray, r{v}.sinArray, rmax] = calculateCoeffs(image2, wholeMask, vv(v), vv(v));
    disp(v)
end
%%
n=16;
for v=1:4
    vv=[2,4,10,20];
        image = im2double(data{n}.greenImage);
    wholeMask = data{n}.gelMask;
    bleachMask = data{n}.bleachSpot;
    refMask = wholeMask-bleachMask;
    ref = sum(sum(image.*refMask))/sum(sum(refMask));
    
    image2 = image-ref;
    r{v}.ID = calcInitDist(image2, wholeMask, r{v}.cosArray, r{v}.sinArray);
    disp(v)
end
%%
subplot(2,3,1)
imagesc(image2.*wholeMask);
subplot(2,3,2)
imagesc(r{1}.ID);
subplot(2,3,3)
imagesc(r{2}.ID);
subplot(2,3,4)
imagesc(r{3}.ID);
subplot(2,3,5)
imagesc(r{4}.ID);
% subplot(2,3,6)
% imagesc(data{16}.IDGrn);
%%
subplot(2,1,1)
imagesc(data{16}.cosArrayGrn);
subplot(2,1,2)
imagesc(data{16}.sinArrayGrn);
colorbar
%%


N = 1024;
M = 1344;
r = zeros(N,M);
%     x= 1.58*[-fliplr(1:data{n}.y) 0 1:(N-floor(data{n}.y)-1)];
%     y= 1.58*[-fliplr(1:data{n}.x) 0 1:(M-floor(data{n}.x)-1)];
y= 1.58*[-fliplr(1:data{16}.y) 0 1:(N-floor(data{16}.y)-1)];
x= 1.58*[-fliplr(1:data{16}.x) 0 1:(M-floor(data{16}.x)-1)];

for i = 1:length(y)
    r(i,:) = sqrt(y(i)^2+x.^2);   
end
% r = zeros(size(image));
theta = zeros(size(image));

for i=1:length(y)
    for j = 1:length(x)
        % calculate theta value at each point in the image
        theta(i,j) = atan2(y(i),x(j))+pi;
    end
end
% deal with theta-singularity at origin
theta(size(image,1)/2+1,size(image,2)/2+1) = 0;

%%
plot(data{39}.time,data{39}.recSpotGrn)
figure
plot(data{39}.time,data{39}.recGelGreen/sum(sum(data{39}.gelMask))+data{39}.refGrn)
figure
plot(data{39}.time,data{39}.recSpotRed)
figure
plot(data{39}.time,data{39}.recGelRed)

%%
%% dataset in thesis (1:43)
type = {'ctrl','cct1','ctrl','ctrl','cct1','cct2','ctrl','cct1','cct2',...
    'ctrl','cct1','cct2','ctrl','cct1','cct2','cct1','cct2','ctrl','cct1',...
    'ctrl','cct1','cct2','cct2','ctrl','ctrl','ctrl','cct1','ctrl','cct1',...
    'ctrl','ctrl','cct1','ctrl','cct1','cct2','ctrl','cct2','ctrl','cct1',...
    'cct2','ctrl','cct1','cct2','ctrl','cct1','cct2','ctrl','cct2','ctrl','cct1','cct2'};
conc = [0 10 0 0 10 10 0 10 10 0 10 10 0 10 10 10 10 0 10 0 10 10 10 0 0 0 ...
    10 0 10 0 0 10 0 10 10 0 5 0 10 5 0 10 5 0 30 40 0 40 0 30 40];

%%
indexCtrl = [];
indexCct1 = [];
indexCct2 = [];
for n=1:43%length(folders)
    if strcmp(data{n}.type,'ctrl')
        indexCtrl = [indexCtrl n];
    elseif strcmp(data{n}.type,'cct1')
        indexCct1 = [indexCct1 n];
    elseif strcmp(data{n}.type,'cct2')
        indexCct2 = [indexCct2 n];
    end
end

indexCtrl = intersect(indexCtrl,good);
indexCct1 = intersect(indexCct1,good);
indexCct2 = intersect(indexCct2,good);