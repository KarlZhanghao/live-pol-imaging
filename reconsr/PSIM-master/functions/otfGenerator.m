function  otf_param =  otfGenerator( na, lambda, corr)
%% reference: sim_algorithm\OtfProvider.java
% Create a new OTF from a (very basic) estimate.
% The curvature factor account for deviation of
% real-world OTFs from the theoretical optimum.
% corr=1 yields an ideal OTF, a=0.2 .. 0.4 are more realistic,
% corr has no effect on cutoff.
% @param na Objectives NA
% @param lambda Emission wavelength (nm)
% @param corr curvature factor, a = [0..1]
%%
% parameters
otf_param.n_pxl = 512;
otf_param.na = na; 
otf_param.lambda = lambda;
otf_param.cutoff = 1000 * (na * 2) / lambda; 
otf_param.cycles_micron = otf_param.cutoff / otf_param.n_pxl ;
% calculate otf
otf_param.vals	= zeros(1,otf_param.n_pxl);
for i = 1 : otf_param.n_pxl
    v = (i-1)/ otf_param.n_pxl ;
    otf_ideal = (2/pi)*(acos(v) - v*sqrt(1-v*v));
    otf_param.vals(1,i) = otf_ideal*(corr^v);
end
otf_param.multi_band = false;
end