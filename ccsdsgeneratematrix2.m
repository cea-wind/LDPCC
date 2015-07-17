function [ Gqc ] = ccsdsgeneratematrix2( H,M,RATE )
%creating the generate matrix of LDPC codes in CCSDS document (version 2,published in 2007)
% input 
%   H: check matrix
%	M: a parameter assign in document CCSDS 131.1-O-2
% 	RATE: information rate
switch(RATE)
    case 1/2
        K=2;
    case 2/3
        K=4;
    case 4/5
        K=8;
end

P = H(1:3*M,end-3*M+1:end);
Q = H(1:3*M,1:M*K);
W = mod(inv_bin(P)*Q,2)';
IMK = eye(M*K);
G = [IMK W];

rowNum = 1:M/4:1+M*K-M/4;
Gqc = zeros(size(G));
Gqc(rowNum,:) = G(rowNum,:);
for m = 1:4*K+12
    for n = 1:M/4-1
        Gqc(rowNum+n,M/4*(m-1)+1:M/4*m) = circshift(G(rowNum,M/4*(m-1)+1:M/4*m),n,2);
    end
end

figure
spy(Gqc);

hold on;
gridG = zeros(size(G));
gridG(M:M:(K-1)*M,:) = 1;
gridG(:,M:M:(K+2)*M) = 1;
spy(gridG,'r',4);
hold off;

end

