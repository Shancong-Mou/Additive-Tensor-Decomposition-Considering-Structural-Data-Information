function xx_ = genSpError1( d1, d2, d3, A)
% Generate sparse error tensor of size d3, d2, d1, using rectangular.
% xx is of size (d1 x d2 x d3). 
% xx_ is of size (d1, d2, d3). 
% where d1 is sample mode
% A is the smaple number where the error appear
 n =size(A,2);
xx_ = zeros(d3, d2, d1);
Err = zeros(d3,d2);

for ii = 1:n
    for i = floor(3*d3/4):floor(3*d3/4)+1
        Err(i,floor(d2/4):floor(d2/4)+1) = 0.1;%abs(randn(1)/10+0.1);
    end
  xx_(:,:,A(ii)) = xx_(:,:,A(ii)) + Err;
%   xx_(:,:,A(ii)) = xx_(:,:,A(ii)) / norm(xx_(:,:,A(ii)), 'fro');
end

% image(xx_(:,:,1),'CDataMapping','scaled')
end