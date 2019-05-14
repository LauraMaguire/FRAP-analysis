%% Make non-normalized recovery and reference curves

bleachRecovery = zeros(2,n);
ref = zeros(2,n);
wholeGel = zeros(2,n);

for t=1:n
%     bleachRecovery(1,t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*bleachSpot));
%     bleachRecovery(2,t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*bleachSpot));
    
    bleachRecovery(1,t) = sum(sum(im2double(GreenImages{1,t}).*bleachSpot));
    bleachRecovery(2,t) = sum(sum(im2double(RedImages{1,t}).*bleachSpot));
    
    
%   ref(1,t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*refSpot));
%     ref(2,t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*refSpot));
    
%     wholeGel(1,t) = sum(sum(imadjust(im2double(GreenImages{1,t}),grnScale).*total));
%     wholeGel(2,t) = sum(sum(imadjust(im2double(RedImages{1,t}),redScale).*total));
    
    wholeGel(1,t) = sum(sum(im2double(GreenImages{1,t}).*total));
    wholeGel(2,t) = sum(sum(im2double(RedImages{1,t}).*total));
    
    disp(t);
end

%% Make various normalized curves

%bleachRefRatio = sum(sum(bleachSpot))/sum(sum(refSpot));
bleachWholeRatio = sum(sum(bleachSpot))/sum(sum(total));

%norm1 = (bleachRecovery./ref)/bleachRefRatio;
norm2 = (bleachRecovery./wholeGel)/bleachWholeRatio;

%%
hold on
plot(time, norm2.','o');
% plot(time, bleachRecovery.','o');
% plot(time, bleachRecovery2.','o');
hold off
%%

figure
plot(time, wholeGel.');

%%
try
[fitresult, ~] = FRAPfitBessel(time(2:end), norm2(1,(2:end)));
a2(1) = fitresult.a;
c2(1) = fitresult.c;
tau2(1) = fitresult.tau;
catch
    disp('could not fit NTF2 curve');
end
try
[fitresult, ~] = FRAPfitBessel(time(2:end), norm2(2,(2:end)));
a2(2) = fitresult.a;
c2(2) = fitresult.c;
tau2(2) = fitresult.tau;
catch
    disp('could not fit mCherry curve');
end

clear fitresult

% [fitresult, ~] = FRAPfitBessel(time(2:end), norm1(1,(2:end)));
% a1(1) = fitresult.a;
% c1(1) = fitresult.c;
% tau1(1) = fitresult.tau;

% [fitresult, ~] = FRAPfitBessel(time(2:end), norm1(2,(2:end)));
% a1(2) = fitresult.a;
% c1(2) = fitresult.c;
% tau1(2) = fitresult.tau;

clear fitresult

%%
close all
clear data data2 bleach bleach2 F GreenImages1 GreenImages2 i image imageg...
    imageb imager ListOfImages1 ListOfImages2 mov n RedImages1 RedImages2...
    ref segment1 segment2 t datapre GreenImagespre RedImagespre ListOfImagesPre...
    segmentpre composite finalRed finalGreen
%%
save([sav '/results.mat']);