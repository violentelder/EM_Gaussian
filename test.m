clear 
close all
clc
%-------------��ȡͼ��-------------
image=imread('images/image_0463.jpg');

%----------------������ͼ��Ԥ����---------------
%ת��ɫ�ʿռ䵽L * a * b *�ռ�
cform = makecform('srgb2lab'); 
lab_i = applycform(image,cform);
figure(1)
subplot(1,2,1)
imshow(lab_i)
title('Ԥ������ͼ��')

% ѡ����ʵ�ͨ�����зָ�
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
%��ֵ�˲�
H1 = medfilt3(image_g);

%---------------�趨EM�㷨�ĳ�ֵ-----------
%���þ�����������Ҫ����Ŀ���
cluster_num =2;
%�����ĳ�ֵ
mu = [107, 150];
%����ĳ�ֵ
sigma = [38^2, 27^2];
%����һ��������ʾÿ����ĸ���ֵ������Ϊ���������������Ϊͼ������ص�ĸ���
pw = zeros(cluster_num,size(H1,1)*size(H1,2));
%�������һ�����ֵ(Ȩ��ֵ)
pc = rand(1,cluster_num);
%���ʹ�һ��
pc = pc/sum(pc);
%����������
max_iter = 20;
iter = 1;
%Ԥ���ռ�
M = zeros(max_iter, cluster_num);
S = zeros(max_iter, cluster_num);

while iter <= max_iter
    %----------E-step------------------
    for i = 1:cluster_num
        %����һ����С���ھ����Ӧ������ֵ�ģ�����Ϊͼ��*���������
        MU = repmat(mu(i),size(H1,1)*size(H1,2),1);
        %��˹ģ��        
        temp = 1/sqrt(2*pi*sigma(i))*exp(-(H1(:)-MU).^2/2/sigma(i));
        %��ֹ����0
        temp(temp<0.000001) = 0.000001;
        pw(i,:) = pc(i) * temp;
    end
    %��һ��
    pw = pw./(repmat(sum(pw),cluster_num,1));
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
%����ֵ�뷽��ĵ���������ʾ����
figure(3)
for i = 1:cluster_num
    plot(M(:,i));
    hold on
end
title('��ֵ�仯����');
figure(4)
for i = 1:cluster_num
    plot(S(:,i));
    hold on
end
title('����仯����');

figure(5)
%Ԫ��������
s_image = cell(1:cluster_num);
rgb_label = repmat(label1, [1,1,3]);
for k=1:cluster_num
    color = image;
    color(rgb_label ~= k-1) = 0;
    s_image{k} = color;
end
subplot(1,2,1)
imshow(s_image{1})
title('��1�е�Ŀ��')
subplot(1,2,2)
imshow(s_image{2})
title('��2�е�Ŀ��')
