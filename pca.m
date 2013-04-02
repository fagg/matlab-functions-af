%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pca.m - Principle Component Analysis via Eigenvector Decomposition
%
% Ashton Fagg (ashton@fagg.id.au) - March 2013.
% 
% Usage: [Y,V,D,mu] = pca(X,Np)
%
% Input
%   - X:  input data
%   - Np: number of projection components
%
% Output
%   - Y:  projected data
%   - V:  principle components (eigenvectors) - optional
%   - D:  vector of eigenvalues - optional
%   - mu: mean of input data - optional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y, varargout] = pca(X,Np)

% X has: N dimensions, M samples -> Fukunaga, pg 39 eqn 2.146
[N,M] = size(X);

if N == Np
    error('PCA:Npeqdims',...
          'Np is the same dimensionality as X. No point doing PCA.');
end

if Np > N
    error('PCA:Npgtdims',...
          'Np is greater than the dimensionality of X. Can''t do PCA.');
end

if Np < 1
    error('PCA:Nplt1','Np must be at least 1.');
end

mu = mean(X,2);
Xhat = bsxfun(@minus, X, mu); % Remove the mean from X

if N > M
    % See Fukunaga book for "Small sample-size problem", pg 39-40:
    S = (Xhat'*Xhat)/N;       % gives m by m, Fukunaga pg 39 eqn 2.147
    
    % Eigenvals/vecs are:
    [V,D] = eig(S);
    V = Xhat * V;             %Fukunaga, pg 39 eqns 2.147 + 2.148
    V = V * inv(sqrt(N*D));   %Fukunaga, pg 40 eqn 2.149
    
    % Normalise and sort in descending order
    D = (N/M) * diag(D);
    [~, idx] = sort(D, 'descend');

else
    % More samples than dimensions, vanilla PCA.
    S = (Xhat*Xhat')./M;
    [V,D] = eig(S);
    % Sort the eigenvalues in descending order
    D = diag(D);
    [~, idx] = sort(D, 'descend');
end


% Trim to Np dimensions.
V = V(:,idx(1:Np));
D = D(idx(1:Np));
Y = V' * Xhat;

% Set up the returns
if nargout == 2
    varargout(1) = {V};
elseif nargout == 3
    varargout(1) = {V};
    varargout(2) = {D};
elseif nargout == 4
    varargout(1) = {V};
    varargout(2) = {D};
    varargout(3) = {mu};
end
end