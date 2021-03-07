clear 
close all
clc
%-------------��ȡͼ��-------------
image=imread('images/2020-Honda-Civic-Type-R-008-2160.jpg');
%ת��ɫ�ʿռ䵽L * a * b *�ռ�
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