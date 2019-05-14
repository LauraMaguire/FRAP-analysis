% Experiment information file.

% File path information

% Deal with path naming differences between Mac and Windows.  I don't know
% if this will work with Linux.
info.expFolderPC = 'Z:\Microscopy';
info.expFolderMac = '/Volumes/houghgrp/Microscopy';

info.baseSavePathPC = 'Z:\Processed Images';
info.baseSavePathMac = '/Volumes/houghgrp/Processed Images';

if ispc % If this computer is a PC
    slash = '\'; % use backslashes along path
    %expFolder = info.expFolderPC; % set base path to experiment
    baseSavePath = info.baseSavePathPC; % set base save path
else % If this computer is a Mac
    slash = '/'; % use forward slashes along path
    %expFolder = info.expFolderMac; % set base path to experiment
    baseSavePath = info.baseSavePathMac; % set base save path
end

% expName: experiment file name
info.expName = '190301_FRAP/190301_FRAP_15.vsi';
% saveFolderOverwrite: name of folder where all results should be saved, if
% not standard.  Leave blank ('') if using standard naming.  Useful for
% older files.
info.saveFolderOverwrite = '';

% Date information
info.day =01; % day of the month that experiment was run
info.month = 03; % month experiment was run
info.year = 2019; % year experiment was run
info.dailyIndex = 6; % indexes order in which experiments were run on that day
info.date = [num2str(info.year) '-' num2str(info.month) '-' ...
    num2str(info.day) '_' num2str(info.dailyIndex)]; % calculated date field

% Microscope information
info.scopeName = 'Olympus'; % Olympus or Nikon?
info.frames = 30; % number of frames
info.nChannels = 2; % number of channels
info.tScale = 120; % seconds per frame
info.xScale = 1.58; % microns per pixel (1.58 is for Olympus, 4x)
info.grnGain = 0; % gain on the green channel
info.redGain = 3; % gain on the red channel
info.grnExp = 10; % exposure time for the green channel (in ms)
info.redExp = 20; % exposure time for the red channel (in ms)
info.scopeNotes = 'post-5s DAPI bleach 40x, 45s til start of exp (4x)'; % extra notes about the microsope settings

% Gel information
info.protein = 'cct2'; % name of protein anchored to gel ('ctrl' for no protein)
info.conc = 40; % protein concentration in mg/mL
info.geo = '1uL gel'; % geometry of gel
info.linker = '30% 29:1 ac:bis'; % MW of linker
info.gelNotes = '6% ac, 2mM LAP'; % extra notes about the gel

% Solution information
info.greenName = 'NTF2-FITC'; % name of green species
info.greenConc = 20; % concentration of green species (uM)
info.redName = 'mCherry'; % name of red species
info.redConc = 20; % concentration of red species (uM)
info.buffer = 'PTB'; % buffer
info.solNotes = ''; % extra notes about solution

% Extra notes
info.extraNotes = '4th chamber cct2 seg. 2, 1 min b/w segments';

% Save all information in a matlab file.
% If saveFolderOverwrite is not a blank string, save in that folder.
% Otherwise, make the standard folder 'year-month-day_dailyIndex'.
if strcmp(info.saveFolderOverwrite, '')
    info.saveFolder = info.date;
else
    info.saveFolder = info.saveFolderOverwrite;
end
fullSavePath = [baseSavePath slash info.saveFolder];
mkdir(fullSavePath);
path = [fullSavePath slash info.date '--info.mat'];
clear expFolder baseSavePath fullSavePath slash
save(path, 'info');

