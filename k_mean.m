clear
close all
%======�����㷨====
%k��ֵ�㷨
i = imread('images/image_0464.jpg');
figure
subplot(2,3,1)
imshow(i)
title('ԭʼͼ��')

%ת��ɫ�ʿռ䵽L * a * b *�ռ�
cform = makecform('srgb2lab'); 
lab_i = applycform(i,cform);

%ʹ��k��ֵ�����㷨���з���
ab = double(lab_i(:,:,2:3));
nrow = size(ab,1);
ncol = size(ab,2);
ab = reshape(ab, nrow * ncol, 2);
%�ظ��������Σ�����ֲ���Сֵ
ncolors = 2;
%ʹ��k��ֵ�����㷨�õ��Ľ�����б��
[c_idx, c_center] = kmeans(ab, ncolors, 'distance', 'sqeuclidean', 'Replicates', 2);
pixel_labels = reshape(c_idx, nrow, ncol);
subplot(2,3,2)
imshow(pixel_labels,[])
title('ʹ�ô�������ͼ����б��')
%Ԫ��������
s_image = cell(1:3);
rgb_label = repmat(pixel_labels, [1,1,3]);
for k=1:ncolors
    color = i;
    color(rgb_label ~= k) = 0;
    s_image{k} = color;
end
subplot(2,3,3)
imshow(s_image{1})
title('��1�е�Ŀ��')
subplot(2,3,4)
imshow(s_image{2})
title('��2�е�Ŀ��')
