function [x,y] = findCenter(mask)


gelRad = zeros(size(mask,1),1);
for i = 1:length(gelRad)
    gelRad(i) = nnz(mask(i,:));
end

[rmax,y] = max(gelRad);
diameter = mask(y,:);
startIndex = find(diameter~=0, 1, 'first');
x = startIndex+rmax/2;


end