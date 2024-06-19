%% PetzRecoveryKraus(K)    Generates the Petz recovery map of the Kraus operator set K 
%   This function has 1 required argument:
%     K: a cell array of Kraus operators, which must be trace-preserving
% 
%   K_tc = PetzRecoveryKraus(K) is a cell array of Kraus operators, which
%   performs the Petz recovery map (or the transpose channel)
% 
%   requires: none
%   author: Bikun Li (bikunli@uchicago.edu)
%   package: QETLAB
%   last updated: June 16, 2024

function K_tc = PetzRecoveryKraus(K)
% % This function converts the input Kraus operator set K as the 
% % Kraus operator set of the Petz recovery map (transpose channel)
% %  version: 06/16/2024

if isempty(K)
    error(['The input ',inputname(1),' is empty'])
else
    K = K(:).';
end

n_K = numel(K);
[n_C,n_S] = size(K{1}); 
% n_C is the dimension of the Hilbert space
% n_S is the dimension of the logical subspace

% % T.P. condition verification
KK = zeros(n_S,n_S);
for i = 1:n_K
    KK = KK + K{i}' * K{i};
end
if norm(KK - eye(n_S),'fro') > 1e-10
    warning(['The input ',inputname(1),' may not be trace-preserving!'])
end

K_mat = cell2mat(K(:).');
M_QEC_mat = full(K_mat'*K_mat); % the QEC matrix
M_QEC_mat_pinv = pinv(M_QEC_mat); % the pseudo inverse of the QEC matrix
M_QEC_mat_pinv_sqrt = sqrtm_homemade(1/2*(M_QEC_mat_pinv + M_QEC_mat_pinv'));
K_tc = mat2cell(M_QEC_mat_pinv_sqrt * (K_mat'),n_S*ones(1,n_K),n_C);

end

function A = sqrtm_homemade(A)
    [V,D] = eig(A);
    A = V*sqrt(D)*V';
end