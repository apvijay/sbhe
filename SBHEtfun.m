% Scrambled Block Hadamard Ensemble (transpose)
%
% Creates and applies \Phi transpose function to a signal based on this paper 
% FAST COMPRESSIVE IMAGING USING SCRAMBLED BLOCK HADAMARD ENSEMBLE
% Gan et al 2008
%
% y        : input vector Mx1
% 
% rowperm  : unif randomly selected M numbers(rows) from 1 to N 
% colperm  : randomly permuted N numbers(columns)
% B        : block size of Hadamard matrix
%
% x        : output vector Nx1
% 
% Example parameter values:
% rowperm = randperm(N,M);
% colperm = randperm(N);
% B = 32;

function x = SBHEtfun(y, rowperm, colperm, B)

M = numel(rowperm);
N = numel(colperm);

% BxB Hadamard matrix
WB = hadamard(B);

x = zeros(M,1);
% Phit = zeros(N,M);

parfor j = 1:N
    colnum = colperm(j);
    
    % choose column, permute based on rowperm, inner prod wit y
    ztop = B * floor((colnum-1)/B);
    phit = [zeros(1,ztop) WB(rem(colnum-1,B)+1,:) zeros(1,N-B-ztop)];
    phit = phit(rowperm);
    x(j,:) = phit * y / sqrt(B); % to make phi unit norm
%     Phit(j,:) = phit;
end
