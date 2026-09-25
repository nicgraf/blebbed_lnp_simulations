%% Radial Variance Transform (RVT)
% This is a MATLAB version of Anna and Alexeys approach to radial symmetry
% we tested this and the difference to the python version is basically
% machine precision (< 1e-14)
%
% written by Anna Kashkanova and Alexey Shkarin
% translated to Matlab by Andre Gemeinhardt
% -- June 2020
%
% inputs: - frame: image to perform algorithm on, does not need to be
%                  square shaped
%         - radii: list of radii (integer) to take into account, they can 
%                  be continuous or with gaps
% outputs: - VoMs: variance of the means map of the image
%          - MoVs: means of the variances map of the image
%
% typically the radial symmetry map is produced by normalizing VoMs by MoVs
% S = VoMs./MoVs, or by just using VoMs as is.
% "The normalized RVT generally is more robust to strong PSF perturbations 
% (such as stretching or PSF astigmatism) and to parameter variations. 
% However, it has relatively high subpixel bias, typically about 0.1 px."

function [VoMs, MoVs] = RVT(frame,radii) 
    frame = frame - mean(frame(:));
    % making arrays
    rmeans = zeros(size(frame,1),size(frame,2),numel(radii));
    rsqmeans = rmeans;
    frame_sq = frame.^2;

    for i = 1:numel(radii)
        % making the circular kernel
        kernel = ringKernel(radii(i),max(radii));
        
        % using MATLAB convolutions for the calculation
        rmeans(:,:,i) = conv2(frame, kernel, 'same');
        rsqmeans(:,:,i) = conv2(frame_sq, kernel, 'same');
    end
    
    rvars = rsqmeans-rmeans.^2; % calculate the variances
    MoVs = mean(rvars,3); % mean of all the variances
    VoMs = var(rmeans,1,3); % variance of all the means
    % note the second argument of the variance being a 1. This is the
    % weight of the normalization. A weight of 0 normalizes by N-1 and a 
    % weight of 1 normalizes by N. Python uses N as well.
end