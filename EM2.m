clear all;close all;clc;
%-------------��ȡͼ��-------------
%img = double(rgb2gray(imread('image_0467.jpg')));
image=imread('image_0467.jpg');
imageg1 = image(:,:,2);%ֻ������ɫͨ��
imageg=double(imageg1);
%----------------������ͼ��Ԥ����---------------
H=fspecial('motion',20,45);
MotionBlur = imfilter(imageg,H,'replicate');
imshow('MotionBlur');
%---------------�趨EM�㷨�ĳ�ֵ-----------
cluster_num =2 ;%���þ�����������Ҫ����Ŀ���
mu = (1:cluster_num)./(cluster_num + 1) .* max(max(H1));%�����ĳ�ֵ
sigma = 20*ones(1,cluster_num);%����ĳ�ֵ
pw = zeros(cluster_num,size(H1,1)*size(H1,2));%����һ��������ʾÿ����ĸ���ֵ������Ϊ���������������Ϊͼ������ص�ĸ���
pc = rand(1,cluster_num);%�������һ�����ֵ(Ȩ��ֵ)  �Ƿ���Բ���������ɣ����Կ�Ч��
pc = pc/sum(pc);%������ʹ�һ��
max_iter = 10;%�Ե�����������Ϊֹͣ������  (�����������ȷ������ʲôӰ��)
iter = 1;
while iter <= max_iter
    %----------E-step------------------
    for i = 1:cluster_num
        MU = repmat(mu(i),size(H1,1)*size(H1,2),1);%ȥmu�еĵ�һ��ֵ������һ��һά����������Ϊ���ظ���
        %��˹ģ��        
        temp = 1/sqrt(2*pi*sigma(i))*exp(-(H1(:)-MU).^2/2/sigma(i));
        temp(temp<0.000001) = 0.000001;%��ֹ����0
        pw(i,:) = pc(i) * temp;
    end
    pw = pw./(repmat(sum(pw),cluster_num,1));%��һ
    %----------M-step--------------------
    %���²�����
    for i = 1:cluster_num
         pc(i) = mean(pw(i,:));
         mu(i) = pw(i,:)*H1(:)/sum(pw(i,:));
         sigma(i) = pw(i,:)*((H1(:)-mu(i)).^2)/sum(pw(i,:));
    end
    %------------show-result---------------
    [~,label] = max(pw);
    %�Ĵ�С
    label = reshape(label,size(H1));
    imshow(label,[])
    title(['iter = ',num2str(iter)]);
    pause(0.1);
    M(iter,:) = mu;
    S(iter,:) = sigma;
    iter = iter + 1;

end
%����ֵ�뷽��ĵ���������ʾ����
figure
for i = 1:cluster_num
    plot(M(:,i));
    hold on
end
title('��ֵ�仯����');
figure
for i = 1:cluster_num
    plot(S(:,i));
    hold on
end
title('����仯����');