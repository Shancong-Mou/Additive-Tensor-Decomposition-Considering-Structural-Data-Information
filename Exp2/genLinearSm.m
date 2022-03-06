function xx_ = genLinearSm(d2, d3, sig)
% Generate circle bkg of size d3, d2.
% xx is of size (d2 x d3). 
% xx_ is of size (d2, d3). 
xx_ = zeros(d3, d2);
for ii = 1:d2
    for jj =1:d3
        xx_(ii,jj) = (ii-1)+(jj-1);
    end
end
 xx_ = xx_/(max(max(xx_)));
% image(xx_(:,:,1),'CDataMapping','scaled')
end