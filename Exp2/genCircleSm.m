function xx_ = genCircleSm(d2, d3)
% Generate circle bkg of size d3, d2.
% xx is of size (d2 x d3). 
% xx_ is of size (d2, d3). 
xx_ = zeros(d3, d2);
sig = 10;
cent = [floor(d2/2);floor(d3/2)];
for ii = 1:d2
    for jj =1:d3
        xx_(ii,jj) = 1/sqrt(2*pi)/sig*exp(-norm([ii;jj]-cent)^2/2/sig^2);
    end
end
 xx_ = xx_/(max(max(xx_)));
% image(xx_(:,:,1),'CDataMapping','scaled')
end