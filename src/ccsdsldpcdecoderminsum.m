function [ iter,decoderData ] = ccsdsldpcdecoderminsum(H,HRowNum_1,HRowNum_2,receiveSignal,MAX_ITER_NUM,alpha)
%CCSDS 中LDPCC的译码算法，（归一化）最小和。
%
% 这一MATLAB函数返回迭代次数和译码结果。
% [ iter,decoderData ] = ccsdsldpcdecoderminsum(H,HRowNum_1,HRowNum_2,receiveSignal,MAX_ITER_NUM,alpha)
% iter ： 迭代次数
% decoderData ： 译码结果
% H：校验矩阵
% receiveSignal : 接收信号
% MAX_ITER_NUM : 限定的最大迭代次数
% alpha：归一化因子
% HRowNum_1,HRowNum_2：更新坐标，可采用如下方式生成
%        [r_mark,~] = find(H~=0);
%         HRowNum_1 = zeros(M,3);
%         for rowH = 1:M
%             HRowNum_1(rowH,:) = find(r_mark==rowH);
%         end
%
%         switch(RATE)
%             case 1/2
%                 HRowNum_2 = zeros(M,6);
%             case 2/3
%                 HRowNum_2 = zeros(M,10);
%             case 3/4
%                 HRowNum_2 = zeros(M,14);
%             case 4/5
%                 HRowNum_2 = zeros(M,18);
%         end
% 
%         for rowH = M+1:3*M
%             HRowNum_2(rowH-M,:) = find(r_mark==rowH);
%         end

[M,N] = size(H);

switch(M/N)
    case 3/5
        numofNz = M*5; colDiv =M/3*cumsum( [0,2,4,2,1,6]);
        vlDiv = M/3*[0,1,2,3,4,5]; rowWeight = [3,6,6];
    case 3/7
        numofNz = M*23/3;colDiv =M/3* cumsum([2*4,2,4,2,1,6]);
        vlDiv = M/3*[2,3,4,5,6,7];rowWeight = [3,10,10];
    case 3/9;
        numofNz = M*31/3;colDiv =M/3* cumsum([4*4,2,4,2,1,6]);
        vlDiv = M/3*[4,5,6,7,8,9];rowWeight = [3,14,14];
    case 3/11
        numofNz = M*13;colDiv =M/3* cumsum([6*4,2,4,2,1,6]);
        vlDiv = M/3*[6,7,8,9,10,11];rowWeight = [3,18,18];
end


vl = receiveSignal;
decoderData = zeros(1,N);

uml = zeros(1,numofNz);
vml1 = reshape(repmat(vl(1:vlDiv(1)),4,1),1,[]);
vml2 = reshape(repmat(vl(vlDiv(1)+1:vlDiv(2)),2,1),1,[]);
vml3 = reshape(repmat(vl(vlDiv(2)+1:vlDiv(3)),4,1),1,[]);
vml4 = reshape(repmat(vl(vlDiv(3)+1:vlDiv(4)),2,1),1,[]);
vml6 = reshape(repmat(vl(vlDiv(5)+1:vlDiv(6)),6,1),1,[]);
vml = [vml1 vml2 vml3 vml4 vl(vlDiv(4)+1:vlDiv(5)) vml6];

for iter=1:MAX_ITER_NUM
    % 行重3   
    vmltemp1 = vml(reshape(HRowNum_1',1,[]));
    vmltemp1 = reshape(vmltemp1,rowWeight(1),[]);
    vmltemp1Mark = ones(size(vmltemp1));
    vmltemp1Mark(vmltemp1<0) = -1;
    vmltemp1Mark_t = repmat(alpha*prod(vmltemp1Mark),rowWeight(1),1);
    vmltemp1Mark = vmltemp1Mark.*vmltemp1Mark_t;
    vmltemp1_min = sort(abs(vmltemp1));
    vmltemp1_min1 = repmat(vmltemp1_min(1,:),rowWeight(1),1);
    vmltemp1_min2 = repmat(vmltemp1_min(2,:),rowWeight(1),1);
    min1_index = find(abs(vmltemp1) == vmltemp1_min1);
    vmltemp1 = vmltemp1_min1;
    vmltemp1(min1_index) = vmltemp1_min2(min1_index);
    vmltemp1 = reshape(vmltemp1.*vmltemp1Mark,1,[]);
    % 行重六
    vmltemp2 = vml(reshape(HRowNum_2',1,[]));
    vmltemp2 = reshape(vmltemp2,rowWeight(2),[]);
    vmltemp2Mark = ones(size(vmltemp2));
    vmltemp2Mark(vmltemp2<0) = -1;
    vmltemp2Mark_t = repmat(alpha*prod(vmltemp2Mark),rowWeight(2),1);
    vmltemp2Mark = vmltemp2Mark.*vmltemp2Mark_t;
    vmltemp2_min = sort(abs(vmltemp2));
    vmltemp2_min1 = repmat(vmltemp2_min(1,:),rowWeight(2),1);
    vmltemp2_min2 = repmat(vmltemp2_min(2,:),rowWeight(2),1);
    min2_index = find(abs(vmltemp2) == vmltemp2_min1);
    vmltemp2 = vmltemp2_min1;
    vmltemp2(min2_index) = vmltemp2_min2(min2_index);
    vmltemp2 = reshape(vmltemp2.*vmltemp2Mark,1,[]);
    uml(reshape(HRowNum_1',1,[])) = vmltemp1;
    uml(reshape(HRowNum_2',1,[])) = vmltemp2;
    %变量节点消息处理
    umltemp1 = reshape(uml(1:colDiv(1)),4,[]);
    umltemp1_sum_vl = sum(umltemp1) + vl(1:vlDiv(1));
    umltemp1 = repmat(umltemp1_sum_vl,4,1) - umltemp1;
    umltemp2 = reshape(uml(colDiv(1)+1:colDiv(2)),2,[]);
    umltemp2_sum_vl = sum(umltemp2) + vl(vlDiv(1)+1:vlDiv(2));
    umltemp2 = repmat(umltemp2_sum_vl,2,1) - umltemp2;
    umltemp3 = reshape(uml(colDiv(2)+1:colDiv(3)),4,[]);
    umltemp3_sum_vl = sum(umltemp3) + vl(vlDiv(2)+1:vlDiv(3));
    umltemp3 = repmat(umltemp3_sum_vl,4,1) - umltemp3;
    umltemp4 = reshape(uml(colDiv(3)+1:colDiv(4)),2,[]);
    umltemp4_sum_vl = sum(umltemp4) + vl(vlDiv(3)+1:vlDiv(4));
    umltemp4 = repmat(umltemp4_sum_vl,2,1) - umltemp4;
    umltemp5 = uml(colDiv(4)+1:colDiv(5));
    umltemp5_sum_vl = umltemp5 + vl(vlDiv(4)+1:vlDiv(5));
    umltemp5 = vl(vlDiv(4)+1:vlDiv(5));
    umltemp6 = reshape(uml(colDiv(5)+1:colDiv(6)),6,[]);
    umltemp6_sum_vl = sum(umltemp6) + vl(vlDiv(5)+1:vlDiv(6));
    umltemp6 = repmat(umltemp6_sum_vl,6,1) - umltemp6;
    
    qn0_1 = [umltemp1_sum_vl umltemp2_sum_vl umltemp3_sum_vl ...
        umltemp4_sum_vl umltemp5_sum_vl umltemp6_sum_vl];
    vml = [reshape(umltemp1,1,[]) reshape(umltemp2,1,[]) reshape(umltemp3,1,[]) ...
        reshape(umltemp4,1,[]) umltemp5 reshape(umltemp6,1,[]) ];
    %译码判决
    decoderData(qn0_1>=0) = 0;
    decoderData(qn0_1<0) = 1;
    if(mod(H*decoderData',2)==0)
        break;
    end
end


