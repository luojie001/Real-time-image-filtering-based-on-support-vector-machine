clc;
clear;
file_path='/home/awen/Awen/training_set/1000';
img_path_list=dir(fullfile(file_path,'*.fit*'));
img_num=length(img_path_list);
for j=1
    image_name=img_path_list(j).name
    image=fitsread(fullfile(file_path,image_name));
    snrf=snrmy(image)
    Source = {10 2.0e11 image 15 3};
    CCD = {1 100 0.99999 100 51 1 10000 0.99*2}; 
    [OutImg,SNRImg] = AstroImgNoise(Source,CCD,0.0001);
    imgname=['/home/awen/Awen/training_set/badsnr',int2str(j+10),'.fits'];
    fitswrite(OutImg,imgname);
    snrimg=snrmy(OutImg)
end
