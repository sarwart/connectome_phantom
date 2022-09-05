function C = my_intersect(A,B)
% Refer to the following URL for this function
%for detailshttps://au.mathworks.com/matlabcentral/answers/51102-ismember-function-too-slow
if ~isempty(A)&&~isempty(B)
 P = zeros(1, max(max(A),max(B)) ) ;
 P(A) = 1;
 C = B(logical(P(B)));
else
  C = [];
end

