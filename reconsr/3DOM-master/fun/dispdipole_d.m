function [dcm,cm, fig] = dispdipole_d(ang,dc,ac,wf_cmin,wf_cmax,img,is_display)
%% 
%%
% wf_cmin = 2000;
% wf_cmax = 17000;
ang = mod(ang,pi);
h = ang/pi;
s = 0.6*ones(size(h));
%% 
wf = dc;
wf = max(wf-wf_cmin, 0); 
wf = min(wf/wf_cmax, 1);
v = wf;
%% 
ld_om_hsv = cat(3,h,s,v);
dcm = hsv2rgb(ld_om_hsv);
% imwrite(dcm,'rho-d.tif');
if nargin < 4
    is_display = false;
end
%
cm = zeros(100, 100);
s = 0.6*ones(size(cm));
%
xx = 0:size(cm,1)-1; yy = xx; c0 = size(cm,1)/2-0.5;
[xx,yy] = meshgrid(xx,yy);
radius = sqrt((xx-c0).^2+(yy-c0).^2);
mask = (radius <= 50) .* (radius>30);  
v = ones(size(cm)).*mask;
%
phy = atan((c0-yy)./(xx-c0));
phy = mod(phy, pi);
h = phy/pi;
%
cm_hsv = cat(3, h, s, v);
cm = hsv2rgb(cm_hsv);
if is_display
    fig = figure;
    dcm = dcm.*(sum(img,3)./max(max(sum(img,3))));
    imshow(dcm),
    title ('azimuthal angle')
    figure;
    imshow(cm)
else
    fig = -1;
end