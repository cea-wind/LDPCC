function [ H ] = ccsdscheckmatrix2( M ,RATE )
%creating the check matrix of LDPC codes in CCSDS document (version 2,published in 2007)
% input 
%	M: a parameter assign in document CCSDS 131.1-O-2
% 	RATE: information rate

load('checkmatrixconstant.mat');

A = zeros(M);
B = eye(M);

L = 0:M-1;
for matrixNum = 1:26
    t_k = theta(matrixNum);
    f_4i_M = floor(4*L/M);
    f_k = fai{matrixNum}(f_4i_M+1,log2(M)-6)';
    col_1 = M/4*(mod((t_k+f_4i_M),4)) + ...
        mod((f_k+L),M/4);
    row_col = col_1+1 + L*M;
    C_temp = zeros(M);
    C_temp(ind2sub([M,M],row_col)) = 1;
    C{matrixNum} = C_temp';
end

H = [A A B A B+C{1};B B A B C{4}+C{3}+C{2};B C{5}+C{6} A C{7}+C{8} B];

switch(RATE)
    case 1/2
        H=H;
        gridCol = 1:4;
    case 2/3
        H_23 = [A A;C{11}+C{10}+C{9} B; B C{14}+C{13}+C{12}];
        H=[H_23 H];
        gridCol = 1:6;
    case 4/5
        H_23 = [A A;C{11}+C{10}+C{9} B; B C{14}+C{13}+C{12}];
        H_45 = [A A A A;C{23}+C{22}+C{21} B C{17}+C{16}+C{15} B;...
            B C{26}+C{25}+C{24} B C{20}+C{19}+C{18}];
        H = [H_45 H_23 H];
        gridCol = 1:10;
end
figure;
spy(H);
hold on;
gridH = zeros(size(H));
gridH([M,2*M],:) = 1;
gridH(:,M*gridCol) = 1;
spy(gridH,'r',4);
hold off;

end

