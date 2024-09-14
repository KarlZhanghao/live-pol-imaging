function sp = simParamCreate(n_band, n_dir, n_phase, size, microns_pxl )
sp.n_dir = n_dir;
sp.n_band  = n_band;
sp.n_phase = n_phase;
sp.pxl  = size;
sp.microns_pixel = microns_pxl;
sp.cycles_micron = 1/(size*microns_pxl);
end