function [rho0, eta0, dc, ac, coeff] = ld3dipole6(pm, alpha, beta)
%%  check input
if ndims(pm) < 3
    pm = reshape(pm,1,1,length(pm));
end
if size(pm,3)~=length(alpha)
    error('data dimension and alpha dimension are unequal !');
end
%% cal ld mat
ldmat = zeros(length(alpha),6);
for kk = 1 : size(ldmat,1)
    ldmat(kk,1) = 1+(sin(beta(kk))).^2;
    ldmat(kk,2) = 2-3*(cos(beta(kk))).^2;
    ldmat(kk,3) = cos(2*alpha(kk)).*(cos(beta(kk))).^2;
    ldmat(kk,4) = sin(2*alpha(kk)).*(cos(beta(kk))).^2;
    ldmat(kk,5) = -2*cos(alpha(kk)).*(sin(2*beta(kk)));
    ldmat(kk,6) = -2*sin(alpha(kk)).*(sin(2*beta(kk)));
end
%
ldmat = ldmat/4;
pmat_inv = pinv(ldmat);
%% cal frequency components
coeff = zeros(size(pm));
for kk = 1 : length(alpha)
    coeff(:,:,1) = coeff(:,:,1) + pmat_inv(1,kk)*pm(:,:,kk);
    coeff(:,:,2) = coeff(:,:,2) + pmat_inv(2,kk)*pm(:,:,kk);
    coeff(:,:,3) = coeff(:,:,3) + pmat_inv(3,kk)*pm(:,:,kk);
    coeff(:,:,4) = coeff(:,:,4) + pmat_inv(4,kk)*pm(:,:,kk);
    coeff(:,:,5) = coeff(:,:,5) + pmat_inv(5,kk)*pm(:,:,kk);
    coeff(:,:,6) = coeff(:,:,6) + pmat_inv(6,kk)*pm(:,:,kk);
end
%% cal dialphae properties
dc = coeff(:,:,1);
ac = sqrt((coeff(:,:,2)).^2+(coeff(:,:,3)).^2+(coeff(:,:,4)).^2+(coeff(:,:,5)).^2+(coeff(:,:,6)).^2);

% 
rho1 = atan(coeff(:,:,4)./coeff(:,:,3))+(coeff(:,:,3)<0).*pi;
rho1(isnan(rho1)) = 0;
rho0 = mod(rho1/2,pi);
%
eta0 =abs((coeff (:,:,6)<0).*pi-atan(sqrt(((coeff(:,:,3)).^2+(coeff(:,:,4)).^2)./((coeff(:,:,5)).^2+(coeff(:,:,6)).^2))));
eta0(isnan(eta0)) = 0;
eta0 = mod(eta0,pi);

end
