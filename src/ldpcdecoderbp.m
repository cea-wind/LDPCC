function [ iter,decoderData ] = ldpcdecoderbp(H,HRowNum,HColNum,receiveSignal,SNR,MAX_ITER_NUM)
%LDPCC decode algorithm , based on the Belief Propagation algorithm
% H: check matrix
% HRowNum,HColNum: index generate from check matrix
% receiveSignal: received soft information from channel
% SNR: signal-to-noise-ration, used to derive metric
% MAT_INTER_NUM: maximum iterations
% HRowNum,HColNum generating by the following codes
% 	[r_mark,c_mark] = find(H~=0);
% 	HColNum = sum(H);
% 	HRowNum = cell(1,size(H,1));
% 	for rowH = 1:size(H,1)
%		 HRowNum{rowH} = find(r_mark==rowH);
% 	end

% The structure of this function is the same as ldpcdecoderminsum


p0 = 1./(1+exp(-2*receiveSignal*SNR));
p1 = 1-p0;
decoderData = zeros(1,length(HColNum));

qnm1 = zeros(1,sum(HColNum));
rmn0 = qnm1;
rmn1 = qnm1;
ColStart = 1;
for L=1:length(HColNum)
    if(p1(L)<eps)
        p1(L) = eps;
        p0(L) = 1-eps;
    end
    if(p0(L)<eps)
        p0(L) = eps;
        p1(L) = 1-eps;
    end
    qnm1(ColStart:ColStart+HColNum(L)-1) = p1(L);
    ColStart = ColStart+HColNum(L);
end

for iter=1:MAX_ITER_NUM
   
    for L=1:length(HRowNum)
        L_col = HRowNum{L};
        qnm1temp = 1-2*qnm1(L_col);temp = qnm1temp;
        for L_col_i = 1:length(L_col)
            v=qnm1temp(L_col_i);
            qnm1temp(L_col_i) = 1;
            qnm1tempprod = prod(qnm1temp);
            if(1-abs(qnm1tempprod)<eps)
                temp(L_col_i) = sign(qnm1tempprod)*(1-2*eps);
            else
                temp(L_col_i) = qnm1tempprod;
            end
            qnm1temp(L_col_i) = v;
        end
        rmn0(L_col)=(1+temp)/2;
        rmn1(L_col)=1-rmn0(L_col);
    end
     
    ColStart = 1;
    qn0_1 = ones(1,length(HColNum));
    for L=1:length(HColNum)
        rmntemp0 = rmn0(ColStart:ColStart+HColNum(L)-1);
        rmntemp1 = rmn1(ColStart:ColStart+HColNum(L)-1);
        temp0 = prod(rmntemp0);
        temp1 = prod(rmntemp1);
        if(abs(temp1)<eps)
            temp1 = sign(1)*eps;
        end
        qn0_1(L) = p0(L)/p1(L)*temp0/temp1;
        if isnan(qn0_1(L))
            disp(qn0_1(L));
            keyboard;
        end
        if(qn0_1(L)==1)
            disp(qn0_1(L));
            keyboard;
        end
        qnm0(ColStart:ColStart+HColNum(L)-1) = p0(L)*temp0./rmntemp0;
        qnm1(ColStart:ColStart+HColNum(L)-1) = p1(L)*temp1./rmntemp1;
        
        k = qnm0(ColStart:ColStart+HColNum(L)-1)+qnm1(ColStart:ColStart+HColNum(L)-1);
        
        qnm1(ColStart:ColStart+HColNum(L)-1) = qnm1(ColStart:ColStart+HColNum(L)-1)./k;
        ColStart = ColStart+HColNum(L);
    end

    decoderData(qn0_1<1) = 1;
    decoderData(qn0_1>=1) = 0;
    if(mod(H*decoderData',2)==0)
        break;
    end
end

