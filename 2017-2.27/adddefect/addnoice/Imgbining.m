function [ Outimg ] = Imgbining( Img,bining )
%Bining ���ͼ��
%   ĿǰΪ���Խ׶Σ�Img��ԭʼͼ��bining�Ƕ��ٵ�Ϊһ�����յ�
Outimg=zeros((size(Img,1)/bining),(size(Img,2)/bining));
for i=1:size(Img,1)/bining
    for j=1:size(Img,1)/bining
        for k=1:bining
            Outimg(i,j)=Outimg(i,j)+Img((i-1)*bining+k,(j-1)*bining+k);
        end
    end
end
end

