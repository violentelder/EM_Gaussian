clear 
close all
clc
%-------------��ȡͼ��-------------
image=imread('images/image_0412.jpg');
%ת��ɫ�ʿռ䵽L * a * b *�ռ�
cform = makecform('srgb2lab'); 
lab_i = applycform(image,cform);
[HIST,~] = histcounts(lab_i(:,:,3),0:1:255, 'Normalization','pdf');
[pks,locs] =findpeaks(HIST,'minpeakdistance',20,'minpeakheight',0.014);
norm(locs)
figure
imhist(lab_i(:,:,1))
figure
imhist(lab_i(:,:,2))
figure
imhist(lab_i(:,:,3))