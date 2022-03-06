%% ATD for example 2
% This code is created by Shancong Mou

function [L,S_s,S_m] = ATDExp2(M,lam,lam1,lam2,lam3,errtol)
%  Parameters:
%       -----------------
%       M      :    Tensor to be decomposed
%       lam    :    ADMM step size
%       lam1   :    Spatial Smoothing Parameter
%       lam2   :    Sparse Parameter for decomposed static hotspot
%       lam3   :    Sparse Parameter for decomposed moving hotspot
%       -----------------
%  Output:
%       -----------------
%       Sp      :    Sparse component
%       Sm      :    Smooth component

%% Default parameters setup
[d3,d2,d1] = size(M);
%% Define matrices
%%
D2 = zeros(d1,d1);
for ii = 2:d1-1
  D2(ii,ii-1) = 1;
  D2(ii,ii) = -2;
  D2(ii,ii+1) = 1;
end
D2(1,1) = -1;
D2(1,2) = 1;
D2(d1,d1-1) = 1;
D2(d1,d1) = -1;


DI1 = zeros(d3-1,d3);
for ii = 1:d3-1
  DI1(ii,ii) = 1;
  DI1(ii,ii+1) = -1;
end

DI2 = zeros(d2-1,d2);
for ii = 1:d2-1
  DI2(ii,ii) = 1;
  DI2(ii,ii+1) = -1;
end

%% tuning parameters
lam = 0.01;
lam1 = 30;
lam2 = 1.9;
lam3 = 2.0;

%% ADMM
d0 = 6;
Z_ = zeros(d3,d2,d1,d0);
U1_ = zeros(d3,d2,d1);
U2_ = zeros(d3,d2,d1);
U3_ = zeros(d3,d2,d1);
U4_ = zeros(d3,d2,d1);
U5_ = zeros(d3,d2,d1);
U6_ = zeros(d3,d2,d1);
A =[1 0 -1 0 0 0; 1 0 0 -1 0 0; 0 1 0 0 -1 0; 1 1 0 0 0 1];
b= zeros(4,d3*d2*d1);
b(4,:)=reshape(M,[],1)';
ATAAT_inv = A'*inv(A*A');
II1 = eye(d3,d3);
II2 = eye(d2,d2);
InvDI1 = inv( 2*lam*lam1* (DI1'*DI1) + II1 );
InvDI2 = inv( 2*lam*lam1* (DI2'*DI2) + II2);
%%
e1 = [];
e2 = [];
while 'true'

Z = Z_;
U1 = U1_;
U2 = U2_;
U3 = U3_;
U4 = U4_;
U5 = U5_;
U6 = U6_;

V1 = Z_(:,:,:,1) - U1;
V2 = Z_(:,:,:,2) - U2;
V3 = Z_(:,:,:,3) - U3;
V4 = Z_(:,:,:,4) - U4;
V5 = Z_(:,:,:,5) - U5;
V6 = Z_(:,:,:,6) - U6;

% l nuc 
[U0,MS,V0] = svd(reshape(V1,d3*d2,d1));
MS_ = max(MS - lam, 0);
S1_ = reshape(U0 * MS_ * V0',d3,d2,d1);
  
[U0,MS,V0] = svd(reshape(V2,d3*d2,d1));
MS_ = max(MS - lam, 0);
S2_ = reshape(U0 * MS_ * V0',d3,d2,d1);
% smooth

for i =1:d1
    S3_(:,:,i) =InvDI1 * (V3(:,:,i));
    
    S4_(:,:,i) = (V4(:,:,i)) * InvDI2;
end
% l_1
S5_ = max(V5-lam*lam2,0)-max(-V5-lam*lam2,0);
S6_ = max(V6-lam*lam3,0)-max(-V6-lam*lam3,0);

%% step 2
X = [reshape(S1_,[],1) reshape(S2_,[],1) reshape(S3_,[],1) reshape(S4_,[],1) reshape(S5_,[],1) reshape(S6_,[],1)];
U = [reshape(U1_,[],1) reshape(U2_,[],1) reshape(U3_,[],1) reshape(U4_,[],1) reshape(U5_,[],1) reshape(U6_,[],1)];

Z0 = ((X+U)'-ATAAT_inv *(A*(X+U)'-b))';
for i =1:6
Z_(:,:,:,i) = reshape(Z0(:,i),d3,d2,d1);
end
U1_ = U1 + S1_ - Z_(:,:,:,1);
U2_ = U2 + S2_ - Z_(:,:,:,2);
U3_ = U3 + S3_ - Z_(:,:,:,3);
U4_ = U4 + S4_ - Z_(:,:,:,4);
U5_ = U5 + S5_ - Z_(:,:,:,5);
U6_ = U6 + S6_ - Z_(:,:,:,6);

Err1 = norm(vec([U1_-U1,U2_-U2,U3_-U3,U4_-U4,U5_-U5,U6_-U6]) ,'fro');
Err2 = norm(vec(Z-Z_),2);
    if (Err1 < errtol)*(Err2 < errtol)
        break
    end
        e1 = [e1;Err1];
e2 = [e2;Err2];

end
    size(e1)
semilogy(e1)
hold on
semilogy(e2)
xlabel('Iteration Number')
ylabel('Residual Error')
legend('Primal residual','Dual residual' )
saveas(gcf,'ResidualPlotExp2.png')
% Decomposed sparse component
L = S1_; % decomposed smooth background
S_s = S5_; % decomposed static hotspot
S_m = S6_; % decomposed moving hotspot
end