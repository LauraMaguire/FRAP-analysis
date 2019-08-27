function rawRecon = calcTimeDist(image,mask,cosArray,sinArray,D,t)

numTerms = size(cosArray,1);
numZeros = size(cosArray,2);

[xCenter,yCenter] = findCenter(mask);

N = size(image,1);
M = size(image,2);
r = zeros(N,M);
y= 1.58*[-fliplr(1:yCenter) 0 1:(N-floor(yCenter)-1)];
x= 1.58*[-fliplr(1:xCenter) 0 1:(M-floor(xCenter)-1)];

for i = 1:length(y)
    r(i,:) = sqrt(y(i)^2+x.^2);   
end

theta = zeros(size(image));
rmax = zeros(length(y),1);

for i = 1:length(y)
    % count nonzero entries in gel mask image
    % largest number of nonzero entries will be the "diameter" in pixels
    rmax(i) = nnz(mask(i,:));
end

for i=1:length(y)
    for j = 1:length(x)
        % calculate theta value at each point in the image
        theta(i,j) = atan2(y(i),x(j))+pi;
    end
end
% deal with theta-singularity at origin
theta(size(image,1)/2+1,size(image,2)/2+1) = 0;

% convert gel "diameter" to "radius" using 1.58 um/pixel
rmax = 1.58*max(rmax)/2;

rawRecon = zeros(size(image));

% fill in several terms so they don't have to get called so often
alpha = zeros(numTerms,numZeros);
cosine = zeros(numTerms,size(image,1),size(image,2));
sine = zeros(numTerms,size(image,1),size(image,2));
jnprimeSq = zeros(numTerms,numZeros);
for n=1:numTerms
    besselOrder = (n-1)-floor(numTerms/2);
    alpha(n,:) = besselzero(besselOrder,numZeros,1)/rmax;
    cosine(n,:,:) = cos(besselOrder.*theta);
    sine(n,:,:) = sin(besselOrder.*theta);
    for a=1:numZeros
        jnprimeSq(n,a) = (0.5*(besselj(besselOrder-1,alpha(n,a)*rmax)...
                -besselj(besselOrder+1,alpha(n,a)*rmax)))^2;  
    end
end

for n=1:numTerms
    besselOrder = (n-1)-floor(numTerms/2);   
    for a=1:numZeros
        jn = (besselj(besselOrder,alpha(n,a)*r));

        termCos = exp(-D.*alpha(n,a).^2.*t).*mask.*squeeze(cosine(n,:,:)).*cosArray(n,a).*jn./jnprimeSq(n,a);
        termSin = exp(-D.*alpha(n,a).^2.*t).*mask.*squeeze(sine(n,:,:)).*sinArray(n,a).*jn./jnprimeSq(n,a);
        rawRecon = rawRecon + termCos + termSin;
        
        disp(['Finished zero number ' num2str(a) ' of ' num2str(numZeros)]);
    end
    disp(['Finished ' num2str(n) ' of ' num2str(numTerms) ' terms.']);
    disp(['Finished Bessel order ' num2str((n-1)-numTerms/2) '.']);
end

end
