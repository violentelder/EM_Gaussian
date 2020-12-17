clear all;close all;clc;
%-------------读取图像-------------
%img = double(rgb2gray(imread('image_0467.jpg')));
image=imread('image_0467.jpg');
imageg1 = image(:,:,2);%只保留绿色通道
imageg=double(imageg1);
%----------------对输入图像预处理---------------
H=fspecial('motion',20,45);
MotionBlur = imfilter(imageg,H,'replicate');
imshow('MotionBlur');
%---------------设定EM算法的初值-----------
cluster_num =2 ;%设置聚类数，既想要分类的块数
mu = (1:cluster_num)./(cluster_num + 1) .* max(max(H1));%期望的初值
sigma = 20*ones(1,cluster_num);%方差的初值
pw = zeros(cluster_num,size(H1,1)*size(H1,2));%构造一个零矩阵表示每个点的概率值，行数为聚类的数量，列数为图像的像素点的个数
pc = rand(1,cluster_num);%随机生成一组概率值(权重值)  是否可以不用随机生成？试试看效果
pc = pc/sum(pc);%将类概率归一化
max_iter = 10;%以迭代次数来作为停止的条件  (迭代次数如何确定？有什么影响)
iter = 1;
while iter <= max_iter
    %----------E-step------------------
    for i = 1:cluster_num
        MU = repmat(mu(i),size(H1,1)*size(H1,2),1);%去mu中的第一个值，构造一个一维向量，行数为像素个数
        %高斯模型        
        temp = 1/sqrt(2*pi*sigma(i))*exp(-(H1(:)-MU).^2/2/sigma(i));
        temp(temp<0.000001) = 0.000001;%防止出现0
        pw(i,:) = pc(i) * temp;
    end
    pw = pw./(repmat(sum(pw),cluster_num,1));%归一
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
    imshow(label,[])
    title(['iter = ',num2str(iter)]);
    pause(0.1);
    M(iter,:) = mu;
    S(iter,:) = sigma;
    iter = iter + 1;

end
%将均值与方差的迭代过程显示出来
figure
for i = 1:cluster_num
    plot(M(:,i));
    hold on
end
title('均值变化过程');
figure
for i = 1:cluster_num
    plot(S(:,i));
    hold on
end
title('方差变化过程');