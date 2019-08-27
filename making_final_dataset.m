% Making the dataset in Thesis Repository/FRAP-analysis-190605.xlsx
% load Google Drive/Hough Lab/Nuclear Pore Team/Laura's Stuff/bound-diffusion-paper/data_final_190605.mat
% The fit results are numFitGrn2 and numFitRed2
% the good indices are [1 2 4:9 11:14 16:28 30:35 37 39:41 43];
% good indices for ctrl are [1,4,7,13,18,20,24,25,26,28,30,31,33,41]
% good indices for cct1 are [2,5,8,11,14,16,19,21,27,32,34,39]
% good indices for cct2 are [6,9,12,17,22,23,35,37,40,43], but last three
% are for 10 mg/mL instead of 20 mg/mL so I left them out
% don't include n > 43 in for loops, those are not good (late additions
% that didn't work well)
% most of the code is in thingIActuallyRan.m (archive folder) but not in a great order