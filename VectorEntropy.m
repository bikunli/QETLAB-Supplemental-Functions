%% VectorEntropy    Compute H = -\sum_i p_i log p_i 
%   This function has 1 required argument:
%     p: a vector for probabilities, which should be normalized so sum(p)=1
% 
%   requires: none
%   author: Bikun Li (bikunli@uchicago.edu)
%   package: QETLAB
%   last updated: June 16, 2024

function [H] = VectorEntropy(p)
% This function compute the value of H_2(p)

% INPUT:
% p can be a scalar (0<=p<=1) of an normalized array
% OUTPUT:
% H: if p is a vector, then H = -\sum_i p_i log p_i
% % if p is a scalar between 0 and 1, then convert p as a vector [p,1-p]

if numel(p) == 1
    p = [p,1-p];
else
    if abs(sum(p(:))-1) > 1e-8
        warning('p is not normalized!')
        p = p/sum(p(:));
    end
    if any(p<0)
        if norm(p(p<0)) > 1e-8
            error('p contains negative entries!');
        end
    end
end

H = 0;
for i = 1:length(p)
    if p(i) > 0
        H = H - p(i)*log2(p(i));
    end
end

end