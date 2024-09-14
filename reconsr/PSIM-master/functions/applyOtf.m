function vec_out = applyOtf( otf, vec_in, kx, ky)
%% reference: sim_algorithm\OtfProvider.java otfToVector
% Multiplies / outputs OTF to a vector. Quite general function,
% some wrappers are provided for conveniece. 
% @param vec  Vector to write / multiply to
% @param band OTF band 
% @param kx OTF center position offset x
% @param ky OTF center position offset y
% @param useAtt if to use attenuation (independent of how {@link #switchAttenuation} is set)
% @param write  if set, vector is overridden instead of multiplied

%%
w = size(vec_in,2);
h = size(vec_in,1);

val = zeros(h,w);

xx = 1: 1: w;
yy = 1: 1: h;
[x,y] = meshgrid(xx,yy);

x(:,1:w/2) = x(:,1:w/2)-1;
x(:,w/2+1:w) = x(:,w/2+1:w) - w - 1;

y(1:h/2,:) = -(y(1:h/2,:)-1);
y(h/2+1:h,:) = h - (y(h/2+1:h,:)-1);

rad = sqrt( (x-kx).^2 + (y-ky).^2 );
cycl = rad * otf.vec_micron;


pos = cycl / otf.cycles_micron + 1;
lPos = floor( pos );
hPos = ceil( pos );
f = pos - lPos;


for i = 1: 1: size(cycl,1)
    for j = 1: 1: size(cycl,2)
        if cycl(i,j) <= otf.cutoff &&  ceil(pos(i,j)) <= otf.n_pxl
            l = 	otf.vals(1,lPos(i,j)) * (1-f(i,j));
            h = 	otf.vals(1,lPos(i,j)) * f(i,j);
            val(i,j) = l + h;
        end
    end
end

vec_out = vec_in.*conj(val);

end