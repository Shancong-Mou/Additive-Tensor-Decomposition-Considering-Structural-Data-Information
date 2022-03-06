%% ATD for example 1
% This code is created by Shancong Mou


function [Sp,Sm] = ATDExp1(M,lam,lam1,lam2,errtol)
%   Parameters:
%       -----------------
%       M      :    Tensor to be decomposed M\in R^d2xd2xd1, where d3 is
%                   the temporal mode     
%       lam    :    ADMM step size
%       lam1   :    Temporal Smoothing Parameter
%       lam2   :    Sparse Parameter
%       -----------------
%  Output:
%       -----------------
%       Sp      :    Sparse component
%       Sm      :    Smooth component

%% Default parameters setup
[d3,d2,d1] = size(M);
%% Define matrices
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
I1 = eye(d1,d1);
II1 = eye(d3,d3);
II2 = eye(d2,d2);
InvDI1 = inv( 2*lam* (DI1'*DI1) + II1 );
InvDI2 = inv( 2*lam* (DI2'*DI2) + II2);
InvDI3 = inv( 2*lam*lam1*(D2'*D2) + I1);
DI1sq = (DI1'*DI1);
DI2sq = (DI2'*DI2);
% e1 = [];
% e2 = [];
%% Use Admm
% initial value 
% lam = 0.01;
% lam1 = 10;
% lam2 = 0.15;

S_bar_ = zeros(d3,d2,d1);
U1_ = zeros(d3,d2,d1);
U2_ = zeros(d3,d2,d1);
U3_ = zeros(d3,d2,d1);
U4_ = zeros(d3,d2,d1);
X1 = reshape(M,d3*d2,d1);
while 'true'
S_bar = S_bar_;
U1 = U1_;
U2 = U2_;
U3 = U3_;
U4 = U4_;
V1 = S_bar - U1;
V2 = S_bar - U2;
V3 = S_bar - U3;
V4 = S_bar - U4;
S1_ = InvDI3*(reshape(V1, d3*d2,d1)');
S1_ = reshape(S1_',d3,d2,d1);

for i =1:d1
    S2_(:,:,i) =InvDI1*(V2(:,:,i) + 2*lam*DI1sq*M(:,:,i));
    S3_(:,:,i) = ( 2*lam*M(:,:,i)* DI2sq +V3(:,:,i)) * InvDI2;
end
S4_ = max(V4-lam*lam2,0)-max(-V4-lam*lam2,0);
S_bar_ = (S1_+ S2_+ S3_+ S4_)/4;
U1_ = U1 + S1_ - S_bar_;
U2_ = U2 + S2_ - S_bar_;
U3_ = U3 + S3_ - S_bar_;
U4_ = U4 + S4_ - S_bar_;
Err1 = norm(reshape(4*(S_bar-S_bar_),d2*d3,d1) ,'fro');
Err2 = norm(vec([U1, U2, U3, U4]-[U1_, U2_, U3_, U4_]),2);
    if (Err1 < errtol)*(Err2<errtol)
        break
    end 
%     e1 = [e1;Err1];
% e2 = [e2;Err2];
end
% size(e1)
% semilogy(e1)
% hold on
% semilogy(e2)
% xlabel('Iteration Number')
% ylabel('Residual Error')
% legend('Primal residual','Dual residual' )
% saveas(gcf,'ResidualPlotExp1.png')
% Decomposed sparse component
Sp = S_bar;
% Decomposed smooth component
Sm = M-Sp;
end