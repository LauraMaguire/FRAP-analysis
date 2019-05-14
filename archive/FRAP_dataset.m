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
%data = cell(length(folders),1);

%% test dataset of high-concentration gels (2-28 and 3-1)
foldersTest = {'19-2-28_1','19-2-28_3','19-2-28_5','19-2-28_7','19-2-28_9',...
    '19-2-28_11','19-2-28_13','19-2-28_15','19-3-1_1','19-3-1_3','19-3-1_5'}.';
data = cell(length(foldersTest),1);
%%
for n=1:length(data)
    r = load(['/Volumes/houghgrp/Processed Images/20' foldersTest{n} '/results.mat']);
    data{n}.grnScale = stretchlim(im2double(r.GreenImages{2}));
    data{n}.redScale = stretchlim(im2double(r.RedImages{2}));
    
    data{n}.greenImage = r.GreenImages{2};
    greenImage = im2double(r.GreenImages{2});
    data{n}.redImage = r.RedImages{2};
    redImage = im2double(r.RedImages{2});
    
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
    
    bleachWholeRatio = sum(sum(bleachSpot))/sum(sum(gelMask));
    data{n}.norm = (bleachRecovery./wholeGel)/bleachWholeRatio;
    
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
    
    data{n}.bProbeq = 1-(data{n}.partCeq(2)/data{n}.partCeq(1));
    data{n}.time = r.time;
    
    disp(['Finished ' num2str(n) ' of ' num2str(length(data)) '.']);
end

clear bleachSpot composite finalGreen finalRed n t gelMask bleachRecovery...
    wholeGel bleachWholeRatio resRef areaPart totRes areaRes totPart greenImageAdj

%%
% for n=14%1:length(data)
%     
%     greenImage = im2double(data{n}.greenImage);
%     greenImage = imadjust(greenImage,data{n}.grnScale);
%     
%     initRed = im2double(data{n}.redImage);
%     initRed = imadjust(initRed,data{n}.redScale);
%     imageb = zeros(size(initRed));
%     
%     image = cat(3,initRed,greenImage,imageb);
%     
%     imshow(image);
%     figure
%     [~, ~, equilMask, ~, ~] = roipoly(image);
%     imshow(equilMask);
%     %data{n}.equilMask = equilMask;
%     close all
%     
%     disp(['Finished ' num2str(n) ' of ' num2str(length(data)) '.']);
% end
% 
% clear bleachSpot composite finalGreen finalRed n t gelMask bleachRecovery...
%     wholeGel bleachWholeRatio resRef areaPart totRes areaRes totPart
%%
for n=1:length(data)
    data{n}.refGrn = sum(sum(im2double(data{n}.greenImage-214).*data{n}.equilMask))...
    /sum(sum(equilMask));
    data{n}.refRed = sum(sum(im2double(data{n}.redImage-219).*data{n}.equilMask))...
    /sum(sum(equilMask));
end
%% dataset in thesis (1:43)
% type = {'ctrl','cct1','ctrl','ctrl','cct1','cct2','ctrl','cct1','cct2',...
%     'ctrl','cct1','cct2','ctrl','cct1','cct2','cct1','cct2','ctrl','cct1',...
%     'ctrl','cct1','cct2','cct2','ctrl','ctrl','ctrl','cct1','ctrl','cct1',...
%     'ctrl','ctrl','cct1','ctrl','cct1','cct2','ctrl','cct2','ctrl','cct1',...
%     'cct2','ctrl','cct1','cct2','ctrl','cct1','cct2','ctrl','cct2','ctrl','cct1','cct2'};
% conc = [0 10 0 0 10 10 0 10 10 0 10 10 0 10 10 10 10 0 10 0 10 10 10 0 0 0 ...
%     10 0 10 0 0 10 0 10 10 0 5 0 10 5 0 10 5 0 30 40 0 40 0 30 40];

%% test dataset
type = {'ctrl','','cct1','cct2','ctrl','cct2','ctrl','cct1','cct2',...
    'ctrl','cct1','cct2'};
conc = [0 30 40 0 40 0 30 40 0 30 40];
%%
for n=1:length(data)
    data{n}.type = type{n};
    data{n}.conc = conc(n);
end
%%
bProb = zeros(length(folders),1);
partC = zeros(length(folders),2);
for n=1:length(folders)
    bProb2(n) = data{n}.bProb2;
    %partC(n,:) = data{n}.partC;
end
%% Collect and sort bound probabilities and partition coefficents
% column 1 is control gels, 2 is cct1, 3 is cct2
% partCR is red partition coefficient, partCG is green
bProbCtrl = [];
partCRCtrl = [];
partCGCtrl = [];
bProbCct1 = [];
partCRCct1 = [];
partCGCct1 = [];
bProbCct2 = [];
partCRCct2 = [];
partCGCct2 = [];

for n=1:length(folders)
    if strcmp(data{n}.type,'ctrl')
        bProbCtrl = [bProbCtrl data{n}.bProb];
        partCRCtrl = [partCRCtrl data{n}.partC(2)];
        partCGCtrl = [partCGCtrl data{n}.partC(1)];
    elseif strcmp(data{n}.type,'cct1')
        bProbCct1 = [bProbCct1 data{n}.bProb];
        partCRCct1 = [partCRCct1 data{n}.partC(2)];
        partCGCct1 = [partCGCct1 data{n}.partC(1)];
    elseif strcmp(data{n}.type,'cct2')
        bProbCct2 = [bProbCct2 data{n}.bProb];
        partCRCct2 = [partCRCct2 data{n}.partC(2)];
        partCGCct2 = [partCGCct2 data{n}.partC(1)];
    end
end

bProb = nan(max([length(bProbCtrl),length(bProbCct1),length(bProbCct2)]),3);
bProb(1:length(bProbCtrl),1) = bProbCtrl;
bProb(1:length(bProbCct1),2) = bProbCct1;
bProb(1:length(bProbCct2),3) = bProbCct2;

partCR = nan(max([length(partCRCtrl),length(partCRCct1),length(partCRCct2)]),3);
partCR(1:length(partCRCtrl),1) = partCRCtrl;
partCR(1:length(partCRCct1),2) = partCRCct1;
partCR(1:length(partCRCct2),3) = partCRCct2;

partCG = nan(max([length(partCGCtrl),length(partCGCct1),length(partCGCct2)]),3);
partCG(1:length(partCGCtrl),1) = partCGCtrl;
partCG(1:length(partCGCct1),2) = partCGCct1;
partCG(1:length(partCGCct2),3) = partCGCct2;

clear bProbCtrl partCRCtrl partCGCtrl bProbCct1 partCRCct1 partCGCct1...
    bProbCct2 partCRCct2 partCGCct2

%% First pass at all fits, stored in 'data' structure.
%     a(n,1) = fitresult.a;
%     c(n,1) = fitresult.c;
%     tau(n,1) = fitresult.tau;
%     ci{n,1} = confint(fitresult);

for n=44:length(folders)
    try
        [fitresult, gof] = FRAPfitBessel(data{n}.time(2:end), ...
            data{n}.norm(1,2:end),n);
        data{n}.fitGrn = fitresult;
        data{n}.gofGrn = gof;
    catch
        disp(['Could not fit NTF2 ' num2str(n)]);
    end
    
    try
        [fitresult, gof] = FRAPfitBessel(data{n}.time(2:end), ...
            data{n}.norm(2,2:end),n);
        data{n}.fitRed = fitresult;
        data{n}.gofRed = gof;
    catch
        disp(['Could not fit mCherry ' num2str(n)]);
    end
end
clear fitresult gof
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
%%
% fitga = nan(length(folders),3);
% fitgc = nan(length(folders),3);
fitgt = nan(length(folders),3);
for n=1:length(folders)
    if ~ismember(n,[5 27 32 35])
        if ismember(n,indexCtrl)
%             fitga(n,1) = data{n}.fitGrn.a;
%             fitgc(n,1) = data{n}.fitGrn.c;
            fitgt(n,1) = data{n}.fitGrn.tau;
       elseif ismember(n,indexCct1)
%             fitga(n,2) = data{n}.fitGrn.a;
%             fitgc(n,2) = data{n}.fitGrn.c;
            fitgt(n,2) = data{n}.fitGrn.tau;
       elseif ismember(n,indexCct2)
%             fitga(n,3) = data{n}.fitGrn.a;
%             fitgc(n,3) = data{n}.fitGrn.c;
            fitgt(n,3) = data{n}.fitGrn.tau;
        end
    end
end
%%
% fitra = nan(length(folders),3);
% fitrc = nan(length(folders),3);
fitrt = nan(length(folders),3);
for n=1:length(folders)
    if ismember(n,indexCtrl)
%         fitra(n,1) = data{n}.fitRed.a;
%         fitrc(n,1) = data{n}.fitRed.c;
        fitrt(n,1) = data{n}.fitRed.tau;
    elseif ismember(n,indexCct1)
%         fitra(n,2) = data{n}.fitRed.a;
%         fitrc(n,2) = data{n}.fitRed.c;
        fitrt(n,2) = data{n}.fitRed.tau;
    elseif ismember(n,indexCct2)
%         fitra(n,3) = data{n}.fitRed.a;
%         fitrc(n,3) = data{n}.fitRed.c;
        fitrt(n,3) = data{n}.fitRed.tau;
    end
end

%% Calculate bleach-spot radius and effective diffusion coefficients
for n=44:length(folders)
    data{n}.w = (1.58)*sqrt(sum(sum(data{n}.bleachSpot))/pi);
    if ~ismember(n,[5 27 32 35 48])
        data{n}.D(1) = data{n}.w^2./data{n}.fitGrn.tau;
    end
    data{n}.D(2) = data{n}.w^2./data{n}.fitRed.tau;
end

%% List and sort green (Dg) and red (Dr) effective diffusion constants
Dg = nan(length(folders),3);
for n=1:length(folders)
    if ~ismember(n,nogood)
        if ismember(n,indexCtrl)
            Dg(n,1) = data{n}.D(1);
        elseif ismember(n,indexCct1)
            Dg(n,2) = data{n}.D(1);
        elseif ismember(n,indexCct2)
            Dg(n,3) = data{n}.D(1);
        end
    end
end
Dr = nan(length(folders),3);
for n=1:length(folders)
    if ismember(n,indexCtrl)
        Dr(n,1) = data{n}.D(2);
    elseif ismember(n,indexCct1)
        Dr(n,2) = data{n}.D(2);
    elseif ismember(n,indexCct2)
        Dr(n,3) = data{n}.D(2);
    end
end
%% First pass at all fits, stored in 'data' structure.
%     a(n,1) = fitresult.a;
%     c(n,1) = fitresult.c;
%     tau(n,1) = fitresult.tau;
%     ci{n,1} = confint(fitresult);

for n=indexCct2
    if ~ismember(n,nogood)

    try
        [fitresult, gof] = FRAPfitBessel(data{n}.time(2:end), ...
            data{n}.norm(1,2:end),n);
        data{n}.fitGrn = fitresult;
        data{n}.gofGrn = gof;
    catch
        disp(['Could not fit NTF2 ' num2str(n)]);
    end
    
    try
        [fitresult, gof] = FRAPfitBessel(data{n}.time(2:end), ...
            data{n}.norm(2,2:end),n);
        data{n}.fitRed = fitresult;
        data{n}.gofRed = gof;
    catch
        disp(['Could not fit mCherry ' num2str(n)]);
    end
    end
end
clear fitresult gof
%%
figure
hold on
for n=indexCtrl
    if n>7
    if ~ismember(n,nogood)
    plot(data{n}.time(2:end), data{n}.norm(1,2:end)/data{n}.norm(1,1),'g-');
    plot(data{n}.time(2:end), data{n}.norm(2,2:end)/data{n}.norm(2,1),'ro');
    end
    end
end
hold off

figure
hold on
for n=indexCct1
    if n>7
    if ~ismember(n,nogood)
    plot(data{n}.time(2:end), data{n}.norm(1,2:end)/data{n}.norm(1,1),'g-');
    plot(data{n}.time(2:end), data{n}.norm(2,2:end)/data{n}.norm(2,1),'ro');
    end
    end
end
hold off

figure
hold on
for n=indexCct2
    if n>7
    if ~ismember(n,nogood)
    plot(data{n}.time(2:end), data{n}.norm(1,2:end)/data{n}.norm(1,1),'g-');
    plot(data{n}.time(2:end), data{n}.norm(2,2:end)/data{n}.norm(2,1),'ro');
    end
    end
end
hold off
%%
ctrlMean=nanmean(fitgt(:,1));
ctrlStdev = nanstd(fitgt(:,1));
cct1Mean=nanmean(fitgt(:,2));
cct1Stdev = nanstd(fitgt(:,2));
cct2Mean=nanmean(fitgt(:,3));
cct2Stdev = nanstd(fitgt(:,3));
%% List and sort goodness-of-fit parameter RMSE
RMSEG = nan(length(folders),3);
for n=1:length(folders)
    if ~ismember(n,[5 27 32 35])
        if ismember(n,indexCtrl)
            RMSEG(n,1) = data{n}.gofGrn.rmse;
        elseif ismember(n,indexCct1)
            RMSEG(n,2) = data{n}.gofGrn.rmse;
        elseif ismember(n,indexCct2)
            RMSEG(n,3) = data{n}.gofGrn.rmse;
        end
    end
end
RMSER = nan(length(folders),3);
for n=1:length(folders)
    if ismember(n,indexCtrl)
        RMSER(n,1) = data{n}.gofRed.rmse;
    elseif ismember(n,indexCct1)
        RMSER(n,2) = data{n}.gofRed.rmse;
    elseif ismember(n,indexCct2)
        RMSER(n,3) = data{n}.gofRed.rmse;
    end
end
%%
semilogy(RMSER,'o');
figure
semilogy(RMSEG,'o');

%% Get upper and lower bounds on tau parameter 
tau_upper = nan(length(folders),3);
for n=1:length(folders)
    if ~ismember(n,nogood)
        ci=confint(data{n}.fitGrn);
        if ismember(n,indexCtrl)
            tau_upper(n,1) = ci(1,3);
        elseif ismember(n,indexCct1)
            tau_upper(n,2) = ci(1,3);
        elseif ismember(n,indexCct2)
            tau_upper(n,3) = ci(1,3);
        end
    end
end
tau_lower = nan(length(folders),3);
for n=1:length(folders)
    if ~ismember(n,nogood)
    ci=confint(data{n}.fitGrn);
    if ismember(n,indexCtrl)
        tau_lower(n,1) = ci(2,3);
    elseif ismember(n,indexCct1)
        tau_lower(n,2) = ci(2,3);
    elseif ismember(n,indexCct2)
        tau_lower(n,3) = ci(2,3);
    end
    end
end
%%
n=3;
semilogy(fitgt(:,n),'o')
hold on
semilogy(tau_upper(:,n),'o')
semilogy(tau_lower(:,n),'o')
hold off
%%
hold on
errorbar(1:43,fitgt(:,1),tau_upper(:,1),tau_lower(:,1),'o')
errorbar(1:43,fitgt(:,2),tau_upper(:,2),tau_lower(:,2),'o')
errorbar(1:43,fitgt(:,3),tau_upper(:,3),tau_lower(:,3),'o')
hold off
%%
chiSquared = zeros(43,1);
reducedChiSquared = zeros(43,1);
for n=1:43
    if ~ismember(n,nogood)
tt = data{n}.time;
yy = data{n}.fitGrn.a.*exp(-data{n}.fitGrn.tau./(2.*tt))...
    .*(besseli(0,data{n}.fitGrn.tau./(2.*tt))+besseli(1,data{n}.fitGrn.tau./(2.*tt)))+data{n}.fitGrn.c;
chiSquare = sum((data{n}.norm(1,2:end)-yy(2:end)).^2./yy(2:end));
chiSquared(n,1) = chiSquare;
rcs = chiSquare/(length(tt(2:end))-2-1);
reducedChiSquared(n,1) = rcs;
    end
end
clear chiSquare rcs tt yy
%%
for n=44:length(folders)
    D(n,:) = data{n}.D;
    bProb2(n) = data{n}.bProb2;
    partC2(n,:) = data{n}.partC;
 end
dBound = (D(:,1)-(1-bProb2').*D(:,2))./bProb2';
%% Get upper and lower bounds on tau parameter 
Dboundsort = nan(length(folders),3);
for n=1:length(folders)
    if ~ismember(n,nogood)
        ci=confint(data{n}.fitGrn);
        if ismember(n,indexCtrl)
            Dboundsort(n,1) = dBound(n);
        elseif ismember(n,indexCct1)
            Dboundsort(n,2) = dBound(n);
        elseif ismember(n,indexCct2)
            Dboundsort(n,3) = dBound(n);
        end
    end
end
%%
for n=26:26
    disp(['loading ' num2str(n)]);
    r = load(['/Volumes/houghgrp/Processed Images/20' folders{n} '/results.mat']);
    areaPart = sum(sum(data{n}.gelMask));
    areaRes = sum(sum(data{n}.resRef));
    
    totPart = sum(sum(im2double(abs(r.GreenImages{1,1}-214)).*data{n}.gelMask));
    totRes = sum(sum(im2double(abs(r.GreenImages{1,1}-219)).*data{n}.resRef));
    data{n}.partC2(1) = (totPart/areaPart)/(totRes/areaRes);
    
    totPart = sum(sum(im2double(abs(r.RedImages{1,1}-214)).*data{n}.gelMask));
    totRes = sum(sum(im2double(abs(r.RedImages{1,1}-219)).*data{n}.resRef));
    data{n}.partC2(2) = (totPart/areaPart)/(totRes/areaRes);
    
    data{n}.bProb2 = 1-(data{n}.partC2(2)/data{n}.partC2(1));
    disp(num2str(n));
end
clear bleachSpot composite finalGreen finalRed n t gelMask bleachRecovery...
    wholeGel bleachWholeRatio resRef areaPart totRes areaRes totPart
%%
%%
for n=2:9
    bk = bfopen(['/Volumes/houghgrp/Microscopy/190226_background/190226_background_0' num2str(n) '.vsi']);
    ListOfImages1 = bk{1,1};
    GreenImages1 = ListOfImages1(1:2:length(ListOfImages1));
    RedImages1 = ListOfImages1(2:2:length(ListOfImages1));
    for t=1:10
        %avggrn(t) = sum(sum(im2double(GreenImages1{1,t})));%/(1024*1344);
        avggrn(t) = sum(sum(GreenImages1{1,t}))/(1024*1344);
        avgred(t) = sum(sum(RedImages1{1,t}))/(1024*1344);
    end
    bkgrndGrn(n) =mean(avggrn);
    bkgrndRed(n) = mean(avgred);
end
