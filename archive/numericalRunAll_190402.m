for n=[39]%length(folders)
    
    [data{n}.x,data{n}.y] = findCenter(data{n}.gelMask);
    
    wholeMask = data{n}.gelMask;
    bleachMask = data{n}.bleachSpot;
    refMask = wholeMask-bleachMask;
    
    image = im2double(data{n}.greenImage);
    data{n}.refGrn = sum(sum(image.*refMask))/sum(sum(refMask));
    image2 = image-data{n}.refGrn;
    
%     [cosArrayGrn, sinArrayGrn, ~] = calculateCoeffs(image2, wholeMask, 2,2);
%     data{n}.cosArrayGrn = cosArrayGrn;
%     data{n}.sinArrayGrn = sinArrayGrn;
%     disp(['Finished green coeffs, ' num2str(n) ' of ' num2str(43)]);
%     
%     image = im2double(data{n}.redImage);
%     data{n}.refRed = sum(sum(image.*refMask))/sum(sum(refMask));
%     image2 = image-data{n}.refRed;
%     
%     [cosArrayRed, sinArrayRed, rmax] = calculateCoeffs(image2, wholeMask, 2,2);
%     data{n}.cosArrayRed = cosArrayRed;
%     data{n}.sinArrayRed = sinArrayRed;
%     data{n}.rmax = rmax;
%     disp(['Finished red coeffs, ' num2str(n) ' of ' num2str(43)]);
%     
%     initialDistribution = calcInitDist(image, wholeMask, cosArrayGrn, sinArrayGrn);
%     data{n}.IDGrn = initialDistribution;
%     disp(['Finished green init. dist., ' num2str(n) ' of 43']);
%     
%     initialDistribution = calcInitDist(image, wholeMask, cosArrayRed, sinArrayRed);
%     data{n}.IDRed = initialDistribution;
%     disp(['Finished red init. dist., ' num2str(n) ' of 43']);
% 
    cosArrayRed=data{n}.cosArrayRed;
    sinArrayRed=data{n}.sinArrayRed;
    cosArrayGrn=data{n}.cosArrayGrn;
    sinArrayGrn=data{n}.sinArrayGrn;
    
%     data{n}.recSpotGrn=simulateData(image,bleachMask,data{n}.cosArrayGrn,...
%         data{n}.sinArrayGrn,data{n}.rmax,data{n}.time,data{n}.numFitGrn.D,data{n}.x,data{n}.y);
%     disp(['Finished green reconstructed spot total, ' num2str(n) ' of 43']);
    
    data{n}.recGelGreen=simulateDataLargeMask(image,wholeMask,cosArrayGrn,...
        sinArrayGrn,data{n}.rmax,data{n}.time,data{n}.numFitGrn.D,data{n}.x,data{n}.y);
    disp(['Finished green reconstructed gel total, ' num2str(n) ' of 43']);
    
    data{n}.recSpotRed=simulateData(image,bleachMask,data{n}.cosArrayRed,...
        data{n}.sinArrayRed,data{n}.rmax,data{n}.time,data{n}.numFitRed.D,data{n}.x,data{n}.y);
    disp(['Finished red reconstructed spot total, ' num2str(n) ' of 43']);
    
    data{n}.recGelRed=simulateDataLargeMask(image,wholeMask,cosArrayRed,...
        sinArrayRed,data{n}.rmax,data{n}.time,data{n}.numFitRed.D,data{n}.x,data{n}.y);
    disp(['Finished red reconstructed gel total, ' num2str(n) ' of 43']);
    
end