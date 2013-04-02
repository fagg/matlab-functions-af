%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pca_recon.m - Reconstruct original data after using Principle Component
% Analysis. See also: pca.m
%
% Ashton Fagg (ashton@fagg.id.au) - April 2013
%
% Xr = pca_recon(Y,V,mu)
%
% Input
%   - Y:   Reduced data
%   - V:   Principle Components (eigenvectors)
%   - mu:  Mean of original data
% Output
%   - Xr:  Reconstructed data
%
% NOTE: The reconstructed data will contain an error, the magnitude of
% which will depend on the number of components used in the projection.
% See also: snr.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Xr = pca_recon(Y,V,mu)
    tmp = V*Y;
    Xr = bsxfun(@plus, tmp, mu);
end