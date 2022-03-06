 % The ATD code is created by Shancong Mou
clc;clear
%import data
load('crack40.mat')
d1=30;d2=40;d3=40;
Sm = genSmTen1( d1, d2, d3, 0.05);
Sp = genSpErrorcrack( d1, d2, d3, [1:d1], a);
M = Sm+Sp;
%% ATD
lam = 0.01;
lam1 = 10;
lam2 = 0.15;
errtol = 1e-6;
tic
[Sp, Sm] = ATDExp1(M,lam,lam1,lam2,errtol);
toc

%% SSD
% The main code for SSD is created by Hao Yan? June 8th, 2017
% If you have any questions, please contact HaoYan@asu.edu
% Paper: Yan, Hao, Kamran Paynabar, and Jianjun Shi. "Anomaly detection in images with smooth background via smooth-sparse decomposition." Technometrics 59.1 (2017): 102-114.
for i =1:d1
dim1 = d3;
dim2 = d2;

spline_degree_bg = 3;
knots_bg1 = dim1; 
B1 = bsplineBasis(dim1,knots_bg1,spline_degree_bg);

knots_bg2 = dim2; 
B2 = bsplineBasis(dim2,knots_bg2,spline_degree_bg);

spline_degree_an_p = 3;
knots_pan1 = dim1;
Bp1 = bsplineBasis(dim1,knots_pan1,spline_degree_an_p);

knots_pan2 = dim2;
Bp2 = bsplineBasis(dim2,knots_pan2,spline_degree_an_p);

bcd_iter = 40;
[Mean_SSD,Anomalies_SSD,lambda_SSD,gamma_SSD,BetaSe] = bsplineSmoothDecompauto(M(:,:,i),{B1,B2},...
                                           {Bp1,Bp2},[],[],bcd_iter);

Anomal(:,:,i) = Anomalies_SSD;
Mean(:,:,i) = Mean_SSD;
end