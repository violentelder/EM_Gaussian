clear
close all
%======聚类算法====
%k均值算法
i = imread('images/image_0464.jpg');
figure
subplot(2,3,1)
imshow(i)
title('原始图像')

%转换色彩空间到L * a * b *空间
cform = makecform('srgb2lab'); 
lab_i = applycform(i,cform);

%使用k均值聚类算法进行分类
ab = double(lab_i(:,:,2:3));
nrow = size(ab,1);
ncol = size(ab,2);
ab = reshape(ab, nrow * ncol, 2);
%重复聚类三次，避免局部最小值
ncolors = 2;
%使用k均值聚类算法得到的结果进行标记
[c_idx, c_center] = kmeans(ab, ncolors, 'distance', 'sqeuclidean', 'Replicates', 2);
pixel_labels = reshape(c_idx, nrow, ncol);
subplot(2,3,2)
imshow(pixel_labels,[])
title('使用簇索引对图像进行标记')
%元胞型数组
s_image = cell(1:3);
rgb_label = repmat(pixel_labels, [1,1,3]);
for k=1:ncolors
    color = i;
    color(rgb_label ~= k) = 0;
    s_image{k} = color;
end
subplot(2,3,3)
imshow(s_image{1})
title('簇1中的目标')
subplot(2,3,4)
imshow(s_image{2})
title('簇2中的目标')
