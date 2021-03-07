clear 
close all
clc
%-------------读取图像-------------
image=imread('images/2020-Honda-Civic-Type-R-008-2160.jpg');
%转换色彩空间到L * a * b *空间
cform = makecform('srgb2lab'); 
lab_i = applycform(image,cform);
[HIST,~] = histcounts(lab_i(:,:,2),0:1:255, 'Normalization','pdf');
[pks,locs] =findpeaks(HIST,'minpeakdistance',20,'minpeakheight',0.014);
norm(locs)
figure
imhist(lab_i(:,:,1))
figure
imhist(lab_i(:,:,2))
figure
imhist(lab_i(:,:,3))