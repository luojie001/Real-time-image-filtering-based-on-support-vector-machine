 function [snr] =snrmy(f)
%求信噪比
x=9;%局部区域所选的大小
[m,n]=size(f);
ph=zeros(m-x,n-x);%用ph来存放 所求的std
parfor k=1:m-x
    for l=1:n-x
    s=f(k:k+x,l:l+x);
     ph(k,l)=std(s(:));
    end
end
phmax=max(max(ph));
phmin=min(min(ph));
snr=10*log(phmax/phmin);
end
