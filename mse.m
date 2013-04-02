%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mse.m - Calculate Mean Square Error (MSE). See also: snr.m
%
% Ashton Fagg (ashton@fagg.id.au) - April 2013
%
% Usage: M = mse(X,Xhat) OR M = mse(E)
%
% Input
%   - X:    Original data
%   - Xhat: Noisy Data
%       OR
%   - E:    Error matrix (from snr.m)
% Output
%   - M:    Mean Square Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = mse(varargin)
    if varargin == 2
        E = varargin{1} - varargin{2};
    end
    if vargin == 1
        E = varargin{1};
    end
    N = size(E,1) * size(E,2);
    M = (1/N) * sum(E.^2);
end

