%% Parameter k
% Variable k should be a matrix of m*2, where m denotes the number of illumination patterns; 
k(1,:) = [131.944,-134.344];
k(2,:) = [-49.144,-181.744];
k(3,:) = [-183.544,-42.011];

%% Parameter phase
% Variable phase should be a vector of 1*m, where m (3 for 2D SIM) denotes the number of illumination patterns;
% phase(i) denotes the referred phase in the ith pattern, the three
% phases in this pattern should be phase(i), phase(i)+2/3*pi and phase(i)+4/3*pi
phase = [-0.349,1.094,-0.288];
%% Directories
read_dir = 'Input\';
saveDir = 'Output\';
calibrationFlag = false; 
if calibrationFlag
    calib1path = 'Calib\calib1.tif'; 
    calib2path = 'Calib\calib2.tif'; 
end

%% whether display the result
displayFlag = true;


%% whether use otf attentuation or not
attFlag = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The PSIM begins
addpath('functions\');
close all;

%% read images
for i = 1: 1: 9
    rawImage(:,:,i) = double(imread([read_dir num2str(i) '.tif'])); 
end

w = size(rawImage,2);
h = size(rawImage,1);

imgs = zeros(size(rawImage));
for i = 1: 1: 9
    imgs(:,:,i) = fadeBorderCos(rawImage(:,:,i),10); 
end

inFFT = zeros(size(rawImage));
for i = 1: 1: 9
    inFFT(:,:,i) = fft2(imgs(:,:,i)); 
end


%% global parameters
showDialog =  false;

nrDirs   = 3;
nrPhases = 3; 
nrBands  = 2;  

lambda = 528;	    
na = 1.4;
pxl_size = 0.080;

otf_corr = 0.31;
wienParam   = 0.05;

otfBeforeShift = true;

findPeak    = true;
refinePhase = false;
	
visualFeedback = 2;
doFastShift = true;

apoB=0.9;
apoF=2;

% SIM
param = simParamCreate(nrBands, nrDirs, nrPhases, w, pxl_size);

% OTF
otf_param = otfGenerator( na, lambda, otf_corr);
otf_param.vec_micron = param.cycles_micron;
otf_raw = oftmatrix( otf_param, w, h);

% Wienar Filter
wiener_filter = wienerGenerator( otf_param, k, w, h);

%% Reconstruction
fullResult = zeros(2*h,2*w);
ld_ft = zeros(h,w,3);

for angIdx = 1 : 3
    kx = k(angIdx,1);
    ky = k(angIdx,2);

    M = [1, 0.5*exp( 1i * (phase(angIdx))), 0.5*exp( -1i * phase(angIdx));
        1, 0.5*exp( 1i * (phase(angIdx)+pi*2/3)), 0.5*exp( -1i * (phase(angIdx)+pi*2/3));
        1, 0.5*exp( 1i * (phase(angIdx)+pi*4/3)), 0.5*exp( -1i * (phase(angIdx)+pi*4/3))];
    invM = inv(M);
    
    separate = zeros(size(inFFT,1),size(inFFT,2),3);
    separate(:,:,1) = invM(1,1) * inFFT(:,:,(angIdx-1)*3+1)+invM(1,2) * inFFT(:,:,(angIdx-1)*3+2)+invM(1,3) * inFFT(:,:,(angIdx-1)*3+3);
    separate(:,:,2) = invM(2,1) * inFFT(:,:,(angIdx-1)*3+1)+invM(2,2) * inFFT(:,:,(angIdx-1)*3+2)+invM(2,3) * inFFT(:,:,(angIdx-1)*3+3);
    separate(:,:,3) = invM(3,1) * inFFT(:,:,(angIdx-1)*3+1)+invM(3,2) * inFFT(:,:,(angIdx-1)*3+2)+invM(3,3) * inFFT(:,:,(angIdx-1)*3+3);
    
    separate_OTF = zeros(size(inFFT,1),size(inFFT,2),3);
    for i=1 : 3  
        separate_OTF(:,:,i) = applyOtf( otf_param, separate(:,:,i), 0, 0);
    end
    
    ld_ft(:,:,angIdx) = separate_OTF(:,:,1);
    
    shifted = zeros(2*w, 2*h,5);
    shifted(:,:,1) = pasteFreq( separate_OTF(:,:,1));
    
    pos = 3;
    neg = 2;
    
    shifted(:,:,pos) = pasteAndFourierShift( separate_OTF(:,:,pos), kx, ky );
    shifted(:,:,neg) = pasteAndFourierShift( separate_OTF(:,:,neg), -kx, -ky );
	
    shifted_mask = zeros(size(shifted));
    shifted_mask(:,:,1) = maskOtf( otf_param, shifted(:,:,1),  0,  0);
    shifted_mask(:,:,pos) = maskOtf( otf_param, shifted(:,:,pos),  kx,  ky);
    shifted_mask(:,:,neg) = maskOtf( otf_param, shifted(:,:,neg),  -kx,  -ky);
    
    thisD = shifted_mask(:,:,1)+shifted_mask(:,:,pos)+shifted_mask(:,:,neg);
    
    
    for i = 1: 1: 3
        fullResult = fullResult + shifted_mask(:,:,i);
    end
end

denom = 1./(wiener_filter+wienParam^2);
fullResult_filtered = fullResult .* denom;

apo = writeApoVector( otf_param, apoB, apoF, 2*h, 2*w);
fullResult_filtered = fullResult_filtered .* apo;


%% pSIM
theta = atan2(k(:,2),k(:,1));
M_ld = 0.5*[1 0.5*exp(2i*theta(1)) 0.5*exp(-2i*theta(1));
    1 0.5*exp(2i*theta(2)) 0.5*exp(-2i*theta(2));
    1 0.5*exp(2i*theta(3)) 0.5*exp(-2i*theta(3))];
invM_ld = inv(M_ld);

p_component(:,:,1) = invM_ld(1,1) * ld_ft(:,:,1) + invM_ld(1,2) * ld_ft(:,:,2) + invM_ld(1,3) * ld_ft(:,:,3);
p_component(:,:,2) = invM_ld(2,1) * ld_ft(:,:,1) + invM_ld(2,2) * ld_ft(:,:,2) + invM_ld(2,3) * ld_ft(:,:,3);
p_component(:,:,4) = invM_ld(3,1) * ld_ft(:,:,1) + invM_ld(3,2) * ld_ft(:,:,2) + invM_ld(3,3) * ld_ft(:,:,3);

psim_f(:,:,1) = zeros(size(fullResult)); psim_f(:,:,3) = fftshift(fullResult);
psim_f(:,:,2) = fftshift(pasteFreq(p_component(:,:,4))); psim_f(:,:,4) = fftshift(pasteFreq(p_component(:,:,2)));

for i = 1: 1: 4
    psim_f(:,:,i) = psim_f(:,:,i).*fftshift(denom);
    psim_f(:,:,i) = psim_f(:,:,i).*fftshift(apo);
end

psim = abs(ifft(ifft(ifft(ifftshift(psim_f),[],1),[],2),[],3));

[sim, psim_om,cm, ouf, ~] = PSIM_display(psim,min(psim(:)),max(psim(:)),false);

%% Display the result 
if displayFlag
    figure;
    imshow(sum(rawImage,3),[]);
    title('WF');
    
    figure;
    imshow(sim,[]);
    title('SIM');
    
    figure;
    imshow(psim_om);
    title('PSIM');
end

%% output result
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

imwrite(uint16(sim), [saveDir 'sim','.tif']);

if calibrationFlag
    imwrite(uint8(psim_om*255), [saveDir 'pSIM_cal', '.png']);
else
    imwrite(uint8(psim_om*255), [saveDir 'pSIM_noncal', '.png']);
end

% Color wheel
imwrite(uint8(cm*255), [saveDir 'cm.png']);

