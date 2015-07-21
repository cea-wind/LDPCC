function [ H ] = ccsdscheckmatrix( M ,RATE )
%creating the check matrix of LDPC codes in CCSDS document (version 1,published in 2006)
% input 
%	M: a parameter assign in document CCSDS 131.1-O-1
% 	RATE: information rate

theta=[1	1	2	3	1	1	2	3	1	2	3	0	2	3	0	2	3	0	1	3	0	1	3	0	1	2];
fai=[1787	1077	1753	697	1523	5	2035	331	1920	130	4	85	551	15	1780	1960	3	145	1019	691	132	42	393	502	201	1064
1502	602	749	1662	1371	9	131	1884	1268	1784	19	1839	81	2031	76	336	529	74	68	186	905	1751	1516	1285	1597	1712
1887	521	590	1775	1738	2032	2047	85	1572	78	26	298	1177	1950	1806	128	1855	129	269	1614	1467	1533	925	1886	2046	1167
1291	301	1353	1405	997	2032	11	1995	623	73	1839	2003	2019	1841	167	1087	2032	388	1385	885	707	1272	7	1534	1965	588];

A = zeros(M);
B = eye(M);

L = 0:M-1;
for matrixNum = 1:26
    t_k = theta(matrixNum);
    f_4i_M = floor(4*L/M);
    f_k = fai(f_4i_M+1,matrixNum)';
    col_1 = M/4*(mod((t_k+f_4i_M),4)) + ...
        mod((f_k+L),M/4);
    row_col = col_1+1 + L*M;
    C_temp = zeros(M);
    C_temp(ind2sub([M,M],row_col)) = 1;
    C{matrixNum} = C_temp';
end

H = [A A A B B+C{1};B+C{8} B+C{7}+C{6} A A B;A B B+C{5} A C{4}+C{3}+C{2}];

switch(RATE)
    case 1/2
        H=H;
        gridCol = 1:4;
    case 2/3
        H_23 = [A A;B C{11}+C{10}+C{9};C{14}+C{13}+C{12} B];
        H=[H_23 H];
        gridCol = 1:6;
    case 4/5
        H_23 = [A A;B C{11}+C{10}+C{9};C{14}+C{13}+C{12} B];
        H_45 = [A A A A;B C{23}+C{22}+C{21} B C{17}+C{16}+C{15};...
            C{26}+C{25}+C{24} B C{20}+C{19}+C{18} B];
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
title('Check Matrix');
hold off;

end

