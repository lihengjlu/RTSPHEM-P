function [outX, outY] = GetSpeRectangles_th20(coarseningFactor)
%obtain side - lengths of rectangles that approximately match the
%permeability distribution provided by SPE10 permeability map
%coarseningFactor^2 data points will be averaged to a single one
% coarseningFactor=300;
global EPS
EPS = eps;
load('color_th20.mat')
Xlength = linspace(0.01, 0.5, 50);
fine = linspace(0.01, 0.5, 200);
[InterpolX, InterpolY] = createLookUp(Xlength, fine);

%temp = Permeability(1:2200, :);
readX = Permeability(1:6900, :)'; 
%temp = Permeability(2201:4400, :);
readY = Permeability(6901:13800, :)';
 
coarseX = zeros(size(readX)/coarseningFactor);
coarseY = zeros(size(readY)/coarseningFactor);

%compute mean value over each batch
for i = 1:(9300 / coarseningFactor)
    for j = 1:(6900 / coarseningFactor)
        temp = readX(((i - 1)*coarseningFactor+1):i*coarseningFactor, ...
            ((j - 1) * coarseningFactor + 1):j*coarseningFactor);
        coarseX(i, j) = mean(temp(:));

        temp = readY(((i - 1)*coarseningFactor+1):i*coarseningFactor, ...
            ((j - 1) * coarseningFactor + 1):j*coarseningFactor);
        coarseY(i, j) = mean(temp(:));
    end
end

%normalize result
coarseX = coarseX / 20000;
coarseY = coarseY / 20000;

% coarseX=log10(coarseX );
% figure()
% imagesc(coarseX );
% shading flat;colorbar;

outX = zeros(size(coarseX));
outY = zeros(size(coarseY));

% map rectangle side-lengths to averaged permeability value by nearest
% neighbor
for i = 1:size(coarseX, 1)
    for j = 1:size(coarseX, 2)
        [~, index] = min((coarseX(i, j)-InterpolX(:)).^2+(coarseY(i, j) - InterpolY(:)).^2);
        temp = meshgrid(fine);
        outX(i, j) = temp(index);
        outY(i, j) = temp(index);
    end
end
%linear projection to regularize
outX = 0.23 + 0.6 * outX;%0.25
outY = 0.23 + 0.6 * outY;

porosity = mean(1-4*outX(:).*outY(:))
end
