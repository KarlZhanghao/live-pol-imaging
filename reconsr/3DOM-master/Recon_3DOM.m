clear all; clc; warning off;close all

%% System parameters
alpha = 0:pi/3:(2*pi-pi/3);
beta = [0.439 0.439 0.439 0.8541 0.8541 0.8541];

%% Read Data
addpath(['.\fun'])
imgdata = bfopen('.\Raw_data.tif');
cordata = bfopen('.\Cor_data.tif');

%% Light intensity correction preprocessing
np = size(cordata{1},1);
mg = zeros(size(cordata{1,1}{1,1}));
x = size(mg,1);
y = size(mg,2);
smg = zeros(x,y,6);
imgstack = zeros(x,y,6);
for pp = 1 : np
    for kk = 1 : 6
        tmp_cor = double(cordata{1}{kk,1});
        suimg(:,:,kk) = max(tmp_cor-100,0);
        smg(:,:,kk)=smg(:,:,kk)+suimg(:,:,kk);
    end
end
smg = smg./np;
wf=mean(smg,3);
for pp = 1 : np
    for kk = 1 : 6
        co(:,:,kk)=wf./double(cordata{1}{kk,1});
        % figure,imshow(co(:,:,kk),[])
        tmp_img = double(imgdata{1}{kk,1});
        tmp_img = max(tmp_img-100,0);
        imgstack(:,:,kk) = tmp_img.*co(:,:,kk); 

    end
end

img0 = double(imgstack(:,:,1:6));

[rho0, eta0, dc, ac, coeff] = ld3dipole6(img0,alpha,beta);
eta0(eta0>(pi/2)) = pi-eta0(eta0>(pi/2));
rho_num = rho0/pi*180;
eta_num = eta0/pi*180;
cmin = 2000;
cmax = 7000;
[dcmrho,cmrho,~] = dispdipole_d(rho0,dc,ac,cmin,cmax,img0,true); 
[dcmeta,cmeta,fig] = dispdipole_eta(eta0,dc,ac,cmin,cmax,img0,true);




