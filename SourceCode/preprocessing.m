clc
clear
close all
%---------读取图片----------%
image = imread('color_th20.tiff');
figure('Name','读取的图片') 
imshow(image) %显示读取的图片
%---------分别提取图片的r，g，b并分别显示rgb和原图----------%
image = imread('color_th20.tiff');
imager = image(:,:,1);
imageg = image(:,:,2);
imageb = image(:,:,3);
figure('Name','图片的RGB图和原图')
subplot(221);
imshow(imager);
title('r')
subplot(222);
imshow(imageg);
title('g')
subplot(223);
imshow(imageb);
title('b')
subplot(224);
imshow(image);
%---------对r，g，b分量加和----------%
figure('Name','图片的RGB图和')
imageSum = imageg+imageb+imager;
imshow(imageSum);%这是rgb加和
[m,n] = size(imageSum);
imageSum2=[6900,9300];
imageSum3=[6900,9300];
parfor i=1:6900
    for j=1:9300
        if imageSum(i,j) >= 255
           imageSum2(i,j) = 0.01;
           imageSum3(i,j) = 0.1;
        else
           imageSum2(i,j) = 0.99;
           imageSum3(i,j) = 1000.0;
        end
    end
end
figure('Name','图片的阈值分割图')
imshow(imageSum2);%这是rgb加和

Porosity=imageSum2;

Permeability=[imageSum3;imageSum3];

save('color_th20','Permeability',"Porosity");