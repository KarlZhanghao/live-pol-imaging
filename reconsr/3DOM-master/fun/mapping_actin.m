function  [rho,eta,mask,disp] = mapping_actin(img_mm,rho_numd,eta_numd,dcmeta,odu,vecZoom,LineWidth)
%% Orientation Mapping
datamat = sum(img_mm,3);
img_ft = fft(img_mm,[],3);
ac_ft = img_ft;
ac_ft(:,:,1) = 0;
acc = ifft(ac_ft,[],3)/size(img_mm,3)*6;
ampl = max(acc,[],3);
th2 = graythresh(ampl./max(ampl(:)));
mask = imbinarize(ampl./max(ampl(:)),th2);


rho_numd1 = rho_numd;
rho_numd1 (mask==0) =0;
rho_numd1(find(rho_numd1 == 0)) = [];
eta_numd1 = eta_numd;
eta_numd1  (mask==0) =0;
eta_numd1(find(eta_numd1 == 0)) = [];
mean_rho = mean(mod(round(180-rho_numd1),180));
eta_numd1(eta_numd1>90) = 180-eta_numd1(eta_numd1>90);
mean_eta = mean(mod(round(eta_numd1),91));
% [x,y] = meshgrid( 1:size(img_mm,2), 1:size(img_mm, 1));
xx = 1 : size(img_mm,2);
yy = 1 : size(img_mm,1);
zz = 1 : size(img_mm,1);
[x,y,z] = meshgrid(xx,yy,zz);
x = x(mask);
y = y(mask);
z = z(mask);

figure,
imshow(datamat,[], 'InitialMagnification','fit'),
title('3D Orientation')
dcmeta0 = zeros(size(dcmeta,1),size(dcmeta,2));
% imshow(dcmeta0, [], 'InitialMagnification','fit'),title('3D Orientation')
hold on;
maxOUF = max(odu(:));
v1 = odu.*cos((180-rho_numd)/180*pi); v1=v1(mask);
u1 = odu.*(sin((180-rho_numd)/180*pi)); u1=u1(mask);
w1 = odu.*(cos((-eta_numd)/180*pi)); w1=w1(mask);
v2 = odu.*cos((-rho_numd)/180*pi); v2=v2(mask);
u2 = odu.*(sin((-rho_numd)/180*pi)); u2=u2(mask);
w2 = odu.*(cos((180-eta_numd)/180*pi)); w2=w2(mask);
dcmeta = abs(dcmeta);
dcmeta1 = dcmeta(:,:,1);dcmeta1 = dcmeta1(mask);
dcmeta2 = dcmeta(:,:,2);dcmeta2 = dcmeta2(mask);
dcmeta3 = dcmeta(:,:,3);dcmeta3 = dcmeta3(mask);
color1 = zeros(size(x));color2 = color1;color3 = color1;
for nn = 1:size(x,1)
    color1(nn) = dcmeta1(nn);
    color2(nn) = dcmeta2(nn);
    color3(nn) = dcmeta3(nn);
end
color1 = gpuArray(color1);
color2 = gpuArray(color2);
color3 = gpuArray(color3);

for kk = 1:size(x,1)
disp = quiver3(x(kk),y(kk),z(kk),v1(kk),u1(kk),w1(kk),vecZoom.*maxOUF, 'color',[color1(kk),color2(kk),color3(kk)], 'LineStyle', '-','LineWidth',LineWidth,'ShowArrowHead','off');%%1.2
disp = quiver3(x(kk),y(kk),z(kk),v2(kk),u2(kk),w2(kk),vecZoom.*maxOUF, 'color',[color1(kk),color2(kk),color3(kk)], 'LineStyle', '-','LineWidth',LineWidth,'ShowArrowHead','off');
end

mean_rho(isnan(mean_rho)) = 0;
mean_eta(isnan(mean_eta)) = 0;
rho = rho_numd;
eta = eta_numd;



end
