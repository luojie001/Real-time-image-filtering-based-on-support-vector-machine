function [OutPutImg,SignalNoiseRatio ] = AstroImgNoise( Source,CCD,Exposuretime )
%AstroImgNoise ��һ����Monte Carlo������ͼ���������ĳ���
%   ����˵���� Source ��һ��cell�����������ԣ�
%               Magnitude: double ��������ǵ�
%               Zeropoint: double  ��ݲ���ȷ���ο���������Բ��ĸ�����
%               Img: matrix ԭʼ��ͼ��
%               SkyBackground: double ���������� ���ǵȣ�
%               TeleA: double ��Զ���ھ���С ���ף�
%%%%%%%%%%%%%%% CCD��һ��cell�����������ԣ�
%               Type: int(1 or 2) 1 ��CCD 2 ��CMOS
%               ReadOutNoise: double CCD�������� (��λ������)
%               ChargeTransferEffciency�� double CCD����ת��Ч�� ���ٷ���
%               DarkCurrent: double CCD������С ����λ���ʣ�
%               Gain: double CCD���� (e/ADU)
%               Bining: int CCD�Ƿ�����bining����ģʽ
%               FullWell: int CCD���С 
%               QuantumEff: double CCD������Ч�� ���ٷ���
%%%%%%%%%%%%%%% ExposureTime: double �ع�ʱ��  ����λ �룩
%I. ������
%Դͼ�������
Mag=cell2mat(Source(1));
Zeropoint=cell2mat(Source(2));
Img=cell2mat(Source(3));
SkyB=cell2mat(Source(4));
TeleA=cell2mat(Source(5));
%CCD������
Type=cell2mat(CCD(1));
RONs=cell2mat(CCD(2));
CTE=cell2mat(CCD(3));
DC=cell2mat(CCD(4));
Gain=cell2mat(CCD(5));
Bining=cell2mat(CCD(6));
FWell=cell2mat(CCD(7));
QE=cell2mat(CCD(8));

%II ������Ͳ����������
if Type == 1
    Gain = Gain*ones(size(Img));
elseif Type == 2
    %CMOS���Ŵ�ЧӦ
    RandAmp=0.01;
    Gain = Gain*(ones(size(Img))+RandAmp*rand((size(Img))));
else
    error('Wrong Type of CCD \n')
end
%III ���ѧ���� ͼ����Ч�Ҷ�ֵת��ΪADU
Nphoton=Zeropoint*10^(-0.4*Mag)*Exposuretime*pi*TeleA^2;
%�������������
Bphoton=Zeropoint*10^(-0.4*SkyB)*Exposuretime*pi*TeleA^2;
%�������طֲ��ĸ����ܶȺ���pdf
Imgpdf=round(Img*Nphoton/(sum(Img(:))));
%��������ͼ�����Ӿ���
Img=Imgpdf*QE+poissrnd(Bphoton/size(Img,1)/size(Img,2),size(Img,1),size(Img,2));
%IV ��������������RON (��Bias�͵�·����,Possion�ֲ�)
RONMatrix=poissrnd(RONs,size(Img,1),size(Img,2));
%�������������� 
DCMatrix=poissrnd(DC*Exposuretime,size(Img,1),size(Img,2));
%��ݵ��ת���ʲ�������ʧ����
CTI=(size(Img,1)+size(Img,2))*(1-CTE).*(1+0.01*rand(size(Img))).*Img;
%�������͵�ɾ�������������ͼ�����
Img=(Img-CTI).*Gain+RONMatrix+DCMatrix;
%���CCD�����������ͼ��
Img((Img>FWell))=FWell;
%IIV ����Ƿ�biningȷ����������С
if (rem(size(Img,1),Bining))||(rem(size(Img,2),Bining))
    error('Can not do bining not correct size \n')
elseif Bining>1
    Img=Imgbining(Img,Bining);
else 
end
%IIIV �������ͼ��
OutPutImg=Img;
%ͼ�������
SignalNoiseRatio=20*log(Nphoton/(...
    +size(Img,1)*size(Img,2)*(RONs^2+DC)+Bphoton+Nphoton)^(0.5));
end

