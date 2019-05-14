function [cosArray, sinArray, rmax] = calculateCoeffs(image, mask, numTerms, numZeros)
tic

[xCenter,yCenter] = findCenter(mask);

N = size(image,1);
M = size(image,2);
r = zeros(N,M);
%     x= 1.58*[-fliplr(1:data{n}.y) 0 1:(N-floor(data{n}.y)-1)];
%     y= 1.58*[-fliplr(1:data{n}.x) 0 1:(M-floor(data{n}.x)-1)];
y= 1.58*[-fliplr(1:yCenter) 0 1:(N-floor(yCenter)-1)];
x= 1.58*[-fliplr(1:xCenter) 0 1:(M-floor(xCenter)-1)];

for i = 1:length(y)
    r(i,:) = sqrt(y(i)^2+x.^2);   
end
% figure
% imshow((r.*mask)/max(max(r)));
% viscircles([xCenter,yCenter],5);

% set up x and y values with origin at center of image
% x=1.58*(-size(image,1)/2:size(image,1)/2-1);
% y=1.58*(-size(image,2)/2:size(image,2)/2-1);

% r = zeros(size(image));
theta = zeros(size(image));
rmax = zeros(length(y),1);

for i = 1:length(y)
    % count nonzero entries in gel mask image
    % largest number of nonzero entries will be the "diameter" in pixels
    rmax(i) = nnz(mask(i,:));
    % calculate radius r at each point in the image
%     r(i,:) = sqrt(x(i)^2+y.^2);
end

for i=1:length(y)
    for j = 1:length(x)
        % calculate theta value at each point in the image
        theta(i,j) = atan2(y(i),x(j))+pi;
    end
end
% deal with theta-singularity at origin
theta(size(image,1)/2+1,size(image,2)/2+1) = 0;
%imshow(theta/max(max(theta)));

% convert gel "diameter" to "radius" using 1.58 um/pixel
rmax = 1.58*max(rmax)/2;
gelArea = pi*(1.58^2)*sum(sum(mask));

cosArray = zeros(numTerms,numZeros);
sinArray = zeros(numTerms,numZeros);

% fill in several terms so they don't have to get called so often
alpha = zeros(numTerms,numZeros);
cosine = zeros(numTerms,size(image,1),size(image,2));
sine = zeros(numTerms,size(image,1),size(image,2));
for n=1:numTerms
    % order nu of the current Bessel function that will be used - integers,
    % centered on zero, total of 'numTerms'
    besselOrder = (n-1)-floor(numTerms/2);
    % list of the first 'numZeros' zeros of all Bessel functions that will
    % be used, normalized to the gel radius rmax
    alpha(n,:) = besselzero(besselOrder,numZeros,1)/rmax;
    % theta-dependence of every integrand, pre-calculated
    cosine(n,:,:) = cos(besselOrder.*theta);
    sine(n,:,:) = sin(besselOrder.*theta);
end

% now do the more time-consuming loop over each term and each zero
for n=1:numTerms
    % calculate the order nu of the current term
    besselOrder = (n-1)-floor(numTerms/2);
    
    % loop over the zeros of that Bessel function
    for a=1:numZeros
        % calculate r-dependence
        jn = (besselj(besselOrder,alpha(n,a).*r));
        jnprimeSq(n,a) = (0.5*real(besselj(besselOrder-1,alpha(n,a)*rmax)...
                -besselj(besselOrder+1,alpha(n,a)*rmax)))^2;  
       
        % create the integrand by multiplying appropriate terms
%         integrandCos = mask.*image.*jn.*squeeze(cosine(n,:,:)).*r;
%         integrandSin = mask.*image.*jn.*squeeze(sine(n,:,:)).*r;
        integrandCos = (1.58)^2*mask.*image.*jn.*squeeze(cosine(n,:,:));
        integrandSin = (1.58)^2*mask.*image.*jn.*squeeze(sine(n,:,:));
        
        intNormalCos = (1.58)^2*mask.*(jn).^2.*(squeeze(cosine(n,:,:))).^2;
        intNormalSin = (1.58)^2*mask.*(jn).^2.*(squeeze(sine(n,:,:))).^2;

        % sum the integrand over the whole image to produce
        % (non-normalized) reconstruction coefficients
%         cosArray(n,a) = (1/gelArea)*sum(sum(integrandCos));
% %         sinArray(n,a) = (1/gelArea)*sum(sum(integrandSin));
%         normCos(n,a) = sum(sum(intNormalCos))/(jnprimeSq(n,a)*gelArea);
%         normSin(n,a) = sum(sum(intNormalSin))/(jnprimeSq(n,a)*gelArea);

        
        normCos(n,a) = sum(sum(intNormalCos))/(jnprimeSq(n,a));
        normSin(n,a) = sum(sum(intNormalSin))/(jnprimeSq(n,a));
        
        cosArray(n,a) = (1/normCos(n,a))*sum(sum(integrandCos));
        sinArray(n,a) = (1/normSin(n,a))*sum(sum(integrandSin));
        if isnan(sinArray(n,a))
            sinArray(n,a) = 0;
        end
        
        %disp(['Finished zero number ' num2str(a) ' of ' num2str(numZeros)]);
    end
    disp(['Finished ' num2str(n) ' of ' num2str(numTerms) ' terms.']);
    disp(['Finished Bessel order ' num2str((n-1)-numTerms/2) '.']);
end

toc
end