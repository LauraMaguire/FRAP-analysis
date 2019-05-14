% after running FRAP_processing_190513, run this script to generate the
% Fourier series fit string and fit the recovery curves

nTerms = 2;

for n=1%:length(data)
    
    [data{n}.x,data{n}.y] = findCenter(data{n}.gelMask);
    
    wholeMask = data{n}.gelMask;
    bleachMask = data{n}.bleachSpot;
    refMask = wholeMask-bleachMask;
    
    image = im2double(data{n}.greenImage);
    %data{n}.refGrn = sum(sum(image.*refMask))/sum(sum(refMask));
    image2 = image-data{n}.refGrn;
    
    [cosArrayGrn, sinArrayGrn, ~] = calculateCoeffs(image2, wholeMask, nTerms,nTerms);
    data{n}.cosArrayGrn = cosArrayGrn;
    data{n}.sinArrayGrn = sinArrayGrn;
    disp(['Finished green coeffs, ' num2str(n) ' of ' num2str(length(data))]);
     
    image = im2double(data{n}.redImage);
    %data{n}.refRed = sum(sum(image.*refMask))/sum(sum(refMask));
    image2 = image-data{n}.refRed;
    
    [cosArrayRed, sinArrayRed, rmax] = calculateCoeffs(image2, wholeMask, nTerms,nTerms);
    data{n}.cosArrayRed = cosArrayRed;
    data{n}.sinArrayRed = sinArrayRed;
    data{n}.rmax = rmax;
    disp(['Finished red coeffs, ' num2str(n) ' of ' num2str(length(data))]);
    
    initialDistribution = calcInitDist(image, wholeMask, cosArrayGrn, sinArrayGrn);
    data{n}.IDGrn = initialDistribution;
    disp(['Finished green init. dist., ' num2str(n) ' of ' num2str(length(data))]);
    
    initialDistribution = calcInitDist(image, wholeMask, cosArrayRed, sinArrayRed);
    data{n}.IDRed = initialDistribution;
    disp(['Finished red init. dist., ' num2str(n) ' of ' num2str(length(data))]);

    cosArrayRed=data{n}.cosArrayRed;
    sinArrayRed=data{n}.sinArrayRed;
    cosArrayGrn=data{n}.cosArrayGrn;
    sinArrayGrn=data{n}.sinArrayGrn;
    
    data{n}.fitStringSpotGrn = generateFitString(image, data{n}.bleachSpot,...
        cosArrayGrn, sinArrayGrn, data{n}.rmax, data{n}.x, data{n}.y);
    
    data{n}.fitStringSpotRed = generateFitString(image, data{n}.bleachSpot,...
        cosArrayRed, sinArrayRed, data{n}.rmax, data{n}.x, data{n}.y);
    
    data{n}.fitStringGelGrn = generateFitString(image, data{n}.gelMask,...
        cosArrayGrn, sinArrayGrn, data{n}.rmax, data{n}.x, data{n}.y);
    
    data{n}.fitStringGelRed = generateFitString(image, data{n}.gelMask,...
        cosArrayRed, sinArrayRed, data{n}.rmax, data{n}.x, data{n}.y);

    
    % insert extra fit parameters to fit to bleach depth and asymptotic
    % recovered value
    data{n}.fGrn = ['(' data{n}.fitStringSpotGrn '/'...
        num2str(sum(sum(data{n}.bleachSpot))) '+' num2str(data{n}.refGrn)...
        ')/(' data{n}.fitStringGelGrn '/' num2str(sum(sum(data{n}.gelMask)))...
        '+' num2str(data{n}.refGrn) ')'];
    
    data{n}.KGrn = ['c1*(' data{n}.fGrn ')+c2'];
    
    data{n}.fRed = ['(' data{n}.fitStringSpotRed '/'...
        num2str(sum(sum(data{n}.bleachSpot))) '+' num2str(data{n}.refRed)...
        ')/(' data{n}.fitStringGelRed '/' num2str(sum(sum(data{n}.gelMask)))...
        '+' num2str(data{n}.refRed) ')'];
    
    data{n}.KRed = ['c1*(' data{n}.fRed ')+c2'];

    [~,data{n}.numFitGrn, ~] = numericalBesselFit(data{n}, data{n}.KGrn,1);

    [~,data{n}.numFitRed, ~] = numericalBesselFit(data{n}, data{n}.KRed,2);

end
clear bleachMask cosArrayGrn cosArrayRed image image2 initialDistribution ...
    n nTerms refMask rmax sinArrayGrn sinArrayRed wholeMask