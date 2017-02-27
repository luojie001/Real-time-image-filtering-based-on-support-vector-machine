clc;
clear;
file_path='/home/awen/Awen/training_set/1000';
img_path_list=dir(fullfile(file_path,'*.fit*'));
img_num=length(img_path_list);
c=ones(1,1000);
%d=ones(1,1000);
e=ones(1,1000);
for j=1:1000
    image_name=img_path_list(j).name;
    image=fitsread(fullfile(file_path,image_name));
    c(1,j)=max(max(image));
    %d(1,j)=sum(sum(image));
    e(1,j)=mean(mean(image));
end