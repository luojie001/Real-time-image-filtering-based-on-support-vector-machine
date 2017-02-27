 % example_ft_sh_phase_screen.m
close all;clear;clc;
file_path='/home/awen/Awen/training_set/1000';
img_path_list=dir(fullfile(file_path,'*.fit*'));
%获取该文件夹中所有的fit格式的图像
img_num=length(img_path_list);
%获取图像总数量
if img_num>0 %满足条件的图像
for jj=1
    image_name=img_path_list(jj).name;%图像名
    image=fitsread(fullfile(file_path,image_name));
    PN=1;
    D = 1; % length of one side of square phase screen [m]
    r0 = 0.1; % coherence diameter [m]
    N = 1024; % number of grid points per side
    L0 = 100; % outer scale [m]
    l0 = 0.01;% inner scale [m]
    % cc=6.88*(D/r0)^(5/3)*2*pi/N;
    delta = D/N; % grid spacing [m]
    % spatial grid
    x = (-N/2 : N/2-1) * delta;
    y = x;
    % generate a random draw of an atmospheric phase screen
    phz= ft_phase_screen(r0, N, delta, L0, l0);
    %screen(:,:,loop)=phz;
    parfor i=1:N
        for j=1:N
            if phz(i,j)>1
                f(i,j)=20*phz(i,j)+300;
            else
                f(i,j)=image(i,j);
            end
        end
    end
    save=['/home/awen/Awen/training_set/mat.fits'];
    fitswrite(f,save);
    %figure, imagesc(f), colorbar
end
end


