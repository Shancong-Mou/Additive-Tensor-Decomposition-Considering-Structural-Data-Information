function xx_ = genSpErrorcrack( d1, d2, d3, A,a)
% Generate sparse error tensor of size d3, d2, d1, using rectangular.
% xx is of size (d1 x d2 x d3). 
% xx_ is of size (d1, d2, d3). 
% where d1 is sample mode
% A is the smaple number where the error appear
 n =size(A,2);
xx_ = zeros(d3, d2, d1);
Err = zeros(d3,d2);


% for ii = 1:n
%     for i = floor(d3/4):floor(d3/2)
%         Err(i,floor(d2/4):floor(d2/4)+1) = abs(randn(1)/10+0.1);
%     end
%   xx_(:,:,A(ii)) = xx_(:,:,A(ii)) + Err;
% %   xx_(:,:,A(ii)) = xx_(:,:,A(ii)) / norm(xx_(:,:,A(ii)), 'fro');
% end
for ii =1 :n
    for i =1:ii
     Err(i,:) = abs(randn(1)/10+0.1)*double(a(i,:));
    end
     xx_(:,:,A(ii)) = xx_(:,:,A(ii)) + Err;
end
% image(xx_(:,:,1),'CDataMapping','scaled')
end