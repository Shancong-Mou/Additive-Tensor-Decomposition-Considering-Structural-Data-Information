function xx_ = genSmTen1( d1, d2, d3, th)
% Generate smooth image of size d2, d3 using gaussian process. 
% xx is of size (d1 x d2 x d3). 
% xx_ is of size (d3, d1, d1). 
% where d1 is sample mode
% for debug: 
% d1 = 30; d2 = 50; th = 0.05;
[x,y] = meshgrid(1:d3,1:d2);
covmat = exp(-th .* pdist2([x(:), y(:)],[x(:), y(:)],'squaredeuclidean'));
xx_ = zeros(d3, d2, d1);
xx = (mvnrnd(zeros(d1, d3 * d2), covmat, d1))';

for ii = 1:d1
xx_(:,:,ii) = reshape(xx(:,ii),  d3, d2);
xx_(:,:,ii) = xx_(:,:,ii) / norm(xx_(:,:,ii), 'fro');
end
% image(xx_(:,:,1),'CDataMapping','scaled')
end