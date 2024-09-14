function wiener_filter = wienerGenerator( otf_param, k, w, h)
%% reference: sim_algorithm\WienerFilter.jave -> updataCache -> addWienerDenominator
% Add OTF^2, for band and direction, to a vector.
% @param d Direction
% @param b Band
% @param useAtt Include attenuation 
wiener_filter = zeros(2*h,2*w);
xx = 1: 1: 2*w;
yy = 1: 1: 2*h;
[x,y] = meshgrid(xx,yy);
x(:,1:w) = x(:,1:w)-1;
x(:,w+1:2*w) = x(:,w+1:2*w)-2*w - 1;

y(1:h,:) = -(y(1:h,:)-1);
y(h+1:2*h,:) = 2*h - (y(h+1:2*h,:)-1);

for d = 1: 1: 3
    for b = 1: 1: 2
        rad1 = sqrt( (x-k(d,1)*(b-1)).^2 + (y-k(d,2)*(b-1)).^2 ) * otf_param.vec_micron;
        rad2 = sqrt( (x+k(d,1)*(b-1)).^2 + (y+k(d,2)*(b-1)).^2 ) * otf_param.vec_micron;
        v_out1 = getOtfVal_m( otf_param, b, rad1);
        v_out2 = getOtfVal_m( otf_param, b, rad2);
        wiener_filter = wiener_filter + v_out1 + v_out2;
    end
end