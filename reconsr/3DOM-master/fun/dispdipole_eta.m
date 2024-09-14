function [dcm,cm, fig] = dispdipole_eta(ang,dc,ac,cmin,cmax,img,is_display)
% h = (pi/2-ang)/(pi/2);
h = (pi-ang)/pi;
s = 0.6*ones(size(h));
%
% v = dc + ac;
% v = v/max(max(v(:))); 
wf = dc;
wf = max(wf-cmin, 0); 
wf = min(wf/cmax, 1)*0.5+0.5; %%将范围从0到1，变成0.5到1
v = wf;
%
ld_om_hsv = cat(3,h,s,v);
dcm = hsv2rgb(ld_om_hsv);

% path = [filepath,'eta_',num2str(eta,'%03d'),'rho_',num2str(rho2,'%03d'),'.tif' ];
% imwrite(dcm,path);

% imwrite(dcm,'eta.tif')
if nargin < 4
    is_display = false;
end
%
cm = zeros(100, 100);
s = 0.6*ones(size(cm));
%
xx = 0:size(cm,1)-1; yy = xx; 
c0 = size(cm,1)/2-0.5;
[xx,yy] = meshgrid(xx,yy);
radius = sqrt((xx-c0).^2+(yy-c0).^2);
mask = (radius <= 50) .* (radius>30);  
% mask = (radius <= 50) .* (radius>0);  

v = ones(size(cm)).*mask;
%
phy = atan((c0-yy)./(xx-c0));
phy = mod(phy, pi);
h = phy/pi;
% phy = mod(phy, pi/2);
% h = phy/(pi/2);


%
cm_hsv = cat(3, h, s, v);
cm = hsv2rgb(cm_hsv);
if is_display
    fig = figure(2);
    dcm = dcm.*(sum(img,3)./max(max(sum(img,3))));
    imshow(dcm)
    title 'polar angle'
%     figure;
%     imshow(cm)
else
    fig = -1;
end