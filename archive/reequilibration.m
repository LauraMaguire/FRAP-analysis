r = load('/Volumes/houghgrp/Processed Images/2019-2-21_11/results.mat');
%%
r2= data{37};
%%
plot(r2.time,r2.norm.')
%%
hold on
plot(r2.time/60,r2.wholeGel(1,:)/r2.wholeGel(1,1),'go')
plot(r2.time/60,r2.wholeGel(2,:)/r2.wholeGel(2,1),'ro')
hold off
xlabel('Time (min)');
ylabel('Entire gel intensity');
legend({'Transport factor','Inert protein'});
%%
finalGreen = im2double(r.GreenImages{1,end});
    finalGreen = imadjust(finalGreen,data{15}.grnScale);
    finalRed = im2double(r.RedImages{1,2});
    finalRed = imadjust(finalRed,data{15}.redScale);
    composite = imfuse(finalGreen, finalRed, 'falsecolor');
    imshow(composite);
    %%
for t=1:length(r.time)
    test(t) = r.GreenImages{1,t}(636,335);
end
%%
plot(test)
%%
