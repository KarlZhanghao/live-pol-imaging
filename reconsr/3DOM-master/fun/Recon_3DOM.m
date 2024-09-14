function [ParaIter,ReconResult] = Recon_3DOM(dataIter)
% setting parameters, lamb, iterations, L
lamb = dataIter.ReconPara(1:2);
iterations = dataIter.ReconPara(3);
L = 0.001;

% read the observed image 
a = dataIter.image;
[M,N,period] = size(a);
imgShape = [M,N];
avgTmp = sum(sum(a))/(M*N);
avg(:) = avgTmp(1,1,:);
modulation = avg/avg(1);
% read psf(point spread function)
psf = dataIter.psf;
psf = psf/sum(sum(psf));
psfShape = size(psf);
offset = ceil(psfShape/2);
% define empty matrix
upImgShape = imgShape;
innerRect = [offset(1),offset(1)+upImgShape(1)-1;...
    offset(2),offset(2)+upImgShape(2)-1];
extendedShape = [ upImgShape(1)+psfShape(1)-1,upImgShape(2)+psfShape(2)-1 ];
% initial calue: g is in spatial domain and 0 initial value
% b is in frequency domain and just not zero in (0,0)
average = sum(sum(sum(a)))/(M*N*period);
b = zeros(imgShape);
b(1,1) = average * sqrt(M*N);
g = zeros( extendedShape(1), extendedShape(2), period );
tkp1 = 1.0;

ParaIter.a = a;
ParaIter.psf = psf;
ParaIter.xk_g = g;
ParaIter.xk_b = b;
ParaIter.ykp1_g = g;
ParaIter.ykp1_b = b;
ParaIter.L = L;
ParaIter.tkp1 = tkp1;
ParaIter.modulation = modulation;
ParaIter.lamb = lamb;
ParaIter.innerRect = innerRect;


% FISTA with backtracking algorithm in [2]
for ii = 1:iterations
   [ParaIter,ReconResult]  = oneIteration(ii, ParaIter);
end
end

function [UpdateParaIter,ReconResult] = oneIteration( ii, ParaIter)
    UpdateParaIter = ParaIter;
    a = ParaIter.a;
    xk_g = ParaIter.xk_g;
    xk_b = ParaIter.xk_b;
    ykp1_g = ParaIter.ykp1_g;
    ykp1_b = ParaIter.ykp1_b;
    psf = ParaIter.psf;
    lamb = ParaIter.lamb;
    L = ParaIter.L;
    tkp1 = ParaIter.tkp1;
    innerRect = ParaIter.innerRect;
    modulation = ParaIter.modulation;
    
    xk1_g = xk_g;
    xk1_b = xk_b;    
    yk_g = ykp1_g;
    yk_b = ykp1_b;
    tk = tkp1;
    % \mu induced by yk_g,yk_b
    forwardProjectionY = forward(yk_g,yk_b,psf,innerRect,modulation); 
    % f(yk)
    maxLikelihoodY = maximumLikelihood(forwardProjectionY, a);
    % grad of f
    [grad1,grad2] = gradient(forwardProjectionY, a, psf, modulation);
    
    for jj = 1:1000
        Ltest = L * power(1.1,jj-1);
        % xtest is p_L(y_k) and step is the operator p_L()
        [xtest_g,xtest_b] = step(Ltest, lamb, yk_g, yk_b, grad1, grad2);
        % new \mu by xtest = p_L(y_k)
        newForwardProjection = forward( xtest_g, xtest_b, psf, innerRect, modulation);
        % f(p_L(y_k)
        newMaxLikelihood = maximumLikelihood(newForwardProjection,a);
        % The first three items of QL()
        quadratic = maxLikelihoodY + quadraticApprox( xtest_g, xtest_b, yk_g, yk_b, grad1, grad2, Ltest);
        if jj > 1
            disp(['Ltest = ', num2str(Ltest)])
            disp(['newMaxLikelihood = ', num2str(newMaxLikelihood)])
            disp(['quadratic = ', num2str(quadratic)])
            disp(['difference = ', num2str(newMaxLikelihood - quadratic)])
        end
        if newMaxLikelihood <= quadratic % inequality condition F(pL(yk))<QL(pL(yk),yk)
            xk_g = xtest_g;
            xk_b = xtest_b;
            L = Ltest;
            maxLikelihoodY = newMaxLikelihood;
            break
        end
        if jj == 10000
            disp('Lipshitz factor still too small after 1000 tries.')
            return
        end
    end
    tkp1 = (1 + sqrt(1 + 4*tk*tk)/2);
    over = (tk -1)/tkp1;
    ykp1_g = xk_g + over * (xk_g- xk1_g);
    ykp1_b = xk_b + over * (xk_b- xk1_b);
    
    UpdateParaIter.xk_g = xk_g;
    UpdateParaIter.xk_b = xk_b;
    UpdateParaIter.ykp1_g = ykp1_g;
    UpdateParaIter.ykp1_b = ykp1_b;
    UpdateParaIter.tkp1 = tkp1;
    UpdateParaIter.L = L;
    
    Mrow = innerRect(1,1):innerRect(1,2);
    Ncol = innerRect(2,1):innerRect(2,2);
    for kk = 1:size(xk_g,3)
        gg(:,:,kk) = xk_g(Mrow,Ncol,kk);
    end
    
    ReconResult.g = gg;
    ReconResult.b = xk_b;

end
function quadraticAdd = quadraticApprox( xtest_g, xtest_b, yk_g, yk_b, grad1, grad2, L)
    delta0 = xtest_g - yk_g;
    delta1 = xtest_b - yk_b;
    quadraticAdd = sum(sum(sum(delta0.*grad1 + L/2*delta0.*delta0)))...
        + sum(sum(sum(delta1.*grad2 + L/2 * delta1 .* delta1)));
end

function [xtest_g,xtest_b] = step(L, lamb, yk_g, yk_b, grad1, grad2)
    c_g = yk_g - grad1/L;
    c_b = yk_b - grad2/L;
    tmp1 = max(c_g - lamb(1)/L, 0);
    tmp2 = min(c_g + lamb(1)/L, 0);
    tmp3 = max(c_b - lamb(2)/L, 0);
    tmp4 = min(c_b + lamb(2)/L, 0);
    xtest_g = max(tmp1 + tmp2, 0);
    xtest_b = tmp3 + tmp4;
end
function [grad1, grad2] = backward(tmp, psf, modulation)
    period = size(tmp,3);
    tmp0 = tmp;
    tmp1 = zeros( size(tmp,1)+size(psf,1)-1,size(tmp,2)+size(psf,2)-1, size(tmp,3));
    tmp2 = zeros(size(tmp,1), size(tmp,2));
    for jj = 1:period
        tmp1(:,:,jj) = modulation(jj) * fftconvolve(tmp0(:,:,jj),psf,'full');
        tmp2 = tmp2 + modulation(jj)*dct2(tmp(:,:,jj));
    end
    grad1 = tmp1;
    grad2 = tmp2;
end
function [grad1,grad2] = gradient(forwardProjectionY, I, psf, modulation);
    tmp = 1 - I./forwardProjectionY;
    [grad1,grad2] = backward(tmp, psf, modulation);
end
function maxLikelihoodY = maximumLikelihood(forwardProjectionY, I)
    temp1 = forwardProjectionY - I .* log(forwardProjectionY);
    maxLikelihoodY = sum(sum(sum(temp1)));
end
function forwardProjectionY = forward(yk_g,yk_b,psf,innerRect,modulation)
    g = yk_g;
    b = yk_b;
    period = size(yk_g,3);
    tmp1 = zeros(size(g));

    Mrow = innerRect(1,1):innerRect(1,2);
    Ncol = innerRect(2,1):innerRect(2,2);
    tmp2 = zeros( size(Mrow,2), size(Ncol,2), period);
    idctB = idct2(b);
    newidctB = zeros([size(idctB),period]);
%     [Mrow,Ncol] = meshgrid(Mrow,Ncol);
    for ii = 1:period
        tmp1(:,:,ii) = modulation(ii) * fftconvolve(g(:,:,ii),psf,'same');
        newidctB(:,:,ii) = modulation(ii) * idctB;
        tmp2(:,:,ii) = tmp1(Mrow,Ncol,ii);
    end
    forwardProjectionY = max( tmp2 + newidctB , 1e-6 );
end
function tmp1 = fftconvolve(g,psf,mode)
    MM = size(g,1) + size(psf,1)-1;
    NN = size(g,2) + size(psf,2)-1;
    G = fft2(g,MM,NN);
    PSF = fft2(psf,MM,NN);
    tmp1 = ifft2(G.*PSF);
    if mode == 'full'
    elseif mode == 'same'
        startpoint1 = ceil((size(psf,1))/2);
        startpoint2 = ceil((size(psf,2))/2);
        tmp2 = tmp1( startpoint1 : startpoint1 + size(g,1)-1, ...
           startpoint2 : startpoint2 + size(g,2) -1 );
       tmp1 = tmp2;
    end
end