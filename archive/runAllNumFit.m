n=39;
data{n}.fitStringSpot = generateFitString(data{n}.greenImage, data{n}.bleachSpot,...
   data{n}.cosArrayGrn, data{n}.sinArrayGrn, data{n}.rmax,...
   data{n}.x,data{n}.y);
data{n}.fitStringGel = generateFitString(data{n}.greenImage, data{n}.gelMask,...
    data{n}.cosArrayGrn, data{n}.sinArrayGrn, data{n}.rmax,...
    data{n}.x,data{n}.y);
data{n}.f = ['(' fitStringSpot '+' num2str(data{n}.refGrn) ')/('...
    fitStringGel '+' num2str(data{n}.refGrn) ')'];
data{n}.K = ['c1*(' f ')+c2'];
[~,data{n}.numFit, data{n}.numGof] = numericalBesselFit(data{n}, K);