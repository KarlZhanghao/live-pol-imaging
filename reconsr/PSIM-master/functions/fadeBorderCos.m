function img_out = fadeBorderCos(img, px)
%% reference: sim_algorithm\SimUtils.java 
% Fades borders (sizes px) of the input to zero.
% Done by multiplying sin^2(x) mapped [0..px] to [pi/2..0].
% Good step before zero-padding data for high-res FFT. 
% @param img Vector to fade borders
% @param px  Size of faded region
%%
w = size(img,2);
h = size(img,1);

img_out = img;

fac = (1/px) *( pi/2);

for y = 1: 1: px
    for x = 1: 1: w
        img_out(y,x) = img_out(y,x) * sin( (y-1) * fac )^2;
    end
end;

for y = h-px+1: 1: h
    for x = 1: 1: w
        img_out(y,x) = img_out(y,x) * sin( (h-y) * fac )^2;
    end
end

for y=1:1:h
    for x = 1: 1: px
        img_out(y,x) = img_out(y,x) * sin( (x-1) * fac )^2;
    end
end

for y=1:1:h
    for x = w-px+1: 1: w
        img_out(y,x) = img_out(y,x) * sin( (w-x) * fac )^2;
    end
end

end