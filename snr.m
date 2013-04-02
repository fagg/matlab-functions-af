%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snr.m - Calculate Signal-to-Noise Ratio between two sets.
%
% Ashton Fagg (ashton@fagg.id.au) - April 2013
%
% Usage: [r,e] = snr(X,Xhat)
%
% Input
%   - X:    Original signal
%   - Xhat: Noisy signal
% Output
%   - r:    SNR in dB
%   - e:    Error matrix (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = snr(X,Xhat)
    e = X - Xhat;
    r = 10*log10(sum(X.^2)/sum(e.^2));
    varargout(1) = {r};
    if (nargout == 2)
        varargout(2) = {e};
    end
end