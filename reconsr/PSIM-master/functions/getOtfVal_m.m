function v_out = getOtfVal_m( otf, band, r_in)

if ~otf.multi_band
    band = 1;
end

v_out = zeros(size(r_in));

pos = r_in / otf.cycles_micron + 1;
lPos = floor( pos );
hPos = ceil( pos );
f = pos - lPos;


for i = 1: 1: size(r_in,1)
    for j = 1: 1: size(r_in,2)
        if r_in(i,j) <= otf.cutoff &&  ceil(pos(i,j)) <= otf.n_pxl
            l = otf.vals(band,lPos(i,j))*(1-f(i,j));
            h = otf.vals(band,lPos(i,j))*f(i,j);
            v_out(i,j) = (l + h)^2;
        end
    end
end

end