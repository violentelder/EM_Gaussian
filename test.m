clear 
close all
clc
%-------------读取图像-------------
image=imread('images/image_0463.jpg');

%----------------对输入图像预处理---------------
%转换色彩空间到L * a * b *空间
cform = makecform('srgb2lab'); 
lab_i = applycform(image,cform);
figure(1)
subplot(1,2,1)
imshow(lab_i)
title('预处理后的图像')

% 选择合适的通道进行分割
aisle = 2;
norm_max = 0;
for i = 1:3
    [HIST,~] = histcounts(lab_i(:,:,i),0:1:255, 'Normalization','pdf');
    [pks,locs] =findpeaks(HIST,'minpeakdistance',20,'minpeakheight',0.014);
    tt = max(locs) - min(locs);
    if(tt > norm_max)
        norm_max = tt;
        aisle = i;
    end
end
image_g = lab_i(:,:,aisle);
image_g = double(image_g);
%中值滤波
H1 = medfilt3(image_g);

%---------------设定EM算法的初值-----------
%设置聚类数，既想要分类的块数
cluster_num =2;
%期望的初值
mu = [107, 150];
%方差的初值
sigma = [38^2, 27^2];
%构造一个零矩阵表示每个点的概率值，行数为聚类的数量，列数为图像的像素点的个数
pw = zeros(cluster_num,size(H1,1)*size(H1,2));
%随机生成一组概率值(权重值)
pc = rand(1,cluster_num);
%概率归一化
pc = pc/sum(pc);
%最大迭代次数
max_iter = 20;
iter = 1;
%预留空间
M = zeros(max_iter, cluster_num);
S = zeros(max_iter, cluster_num);

while iter <= max_iter
    %----------E-step------------------
    for i = 1:cluster_num
        %创建一个大小等于聚类对应期望初值的，长度为图像长*宽的列向量
        MU = repmat(mu(i),size(H1,1)*size(H1,2),1);
        %高斯模型        
        temp = 1/sqrt(2*pi*sigma(i))*exp(-(H1(:)-MU).^2/2/sigma(i));
        %防止出现0
        temp(temp<0.000001) = 0.000001;
        pw(i,:) = pc(i) * temp;
    end
    %归一化
    pw = pw./(repmat(sum(pw),cluster_num,1));
    %----------M-step--------------------
    %更新参数集
    for i = 1:cluster_num
         pc(i) = mean(pw(i,:));
         mu(i) = pw(i,:)*H1(:)/sum(pw(i,:));
         sigma(i) = pw(i,:)*((H1(:)-mu(i)).^2)/sum(pw(i,:));
    end
    %------------show-result---------------
    [~,label] = max(pw);
    %改大小
    label = reshape(label,size(H1));
    figure(1)
    subplot(1,2,2)
    label1 = imbinarize(label, 'adaptive');
    imshow(label1,[])
    title(['iter = ',num2str(iter)]);
    pause(0.1);
    M(iter,:) = mu;
    S(iter,:) = sigma;
    iter = iter + 1;
end
%将均值与方差的迭代过程显示出来
figure(3)
for i = 1:cluster_num
    plot(M(:,i));
    hold on
end
title('均值变化过程');
figure(4)
for i = 1:cluster_num
    plot(S(:,i));
    hold on
end
title('方差变化过程');

figure(5)
%元胞型数组
s_image = cell(1:cluster_num);
rgb_label = repmat(label1, [1,1,3]);
for k=1:cluster_num
    color = image;
    color(rgb_label ~= k-1) = 0;
    s_image{k} = color;
end
subplot(1,2,1)
imshow(s_image{1})
title('簇1中的目标')
subplot(1,2,2)
imshow(s_image{2})
title('簇2中的目标')
