function [ iter,decoderData ] = ldpcdecoderllr(H,HRowNum,HColNum,receiveSignal,SNR,MAX_ITER_NUM)
%LDPCC decode algorithm , based on the log-likelihood ratio algorithm
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

[~,N] = size(H);    

vl = 2*receiveSignal*SNR;
decoderData = zeros(1,N);

uml = zeros(1,sum(HColNum));
vml = uml;
ColStart = 1;
for L=1:length(HColNum)
    vml(ColStart:ColStart+HColNum(L)-1) = vl(L);
    ColStart = ColStart+HColNum(L);
end

for iter=1:MAX_ITER_NUM
 
    for L=1:length(HRowNum)
        L_col = HRowNum{L};
        vmltemp = vml(L_col);
        vmlMark = ones(size(vmltemp));
        vmlMark(vmltemp<0) = -1;
        vmlMark = prod(vmlMark);
        vmltemp1 = (exp(abs(vmltemp))+1)./(exp(abs(vmltemp))-1+eps);
        faitemp = sum(log(vmltemp1+eps));
        for L_col_i = 1:length(L_col)
            vmltemp2 = (exp(abs(vmltemp(L_col_i)))+1)/...
                (exp(abs(vmltemp(L_col_i)))-1+eps);
            vmltemp2 = log(vmltemp2+eps);
            vmltemp3 = faitemp-vmltemp2;
            vmltemp3 = log((exp(vmltemp3)+1)/(exp(vmltemp3)-1+eps)+eps);
            if vmltemp(L_col_i)<0
            	vmltemp(L_col_i) = -vmlMark*vmltemp3;
            else
            	vmltemp(L_col_i) = vmlMark*vmltemp3;
            end
        end
        uml(L_col) = vmltemp;
    end
 
    ColStart = 1;
    qn0_1 = ones(1,N);
    for L=1:length(HColNum)
        umltemp = uml(ColStart:ColStart+HColNum(L)-1);
        temp = sum(umltemp);
        qn0_1(L) = temp + vl(L); 
        umltemp = temp - umltemp;
        vml(ColStart:ColStart+HColNum(L)-1) = umltemp + vl(L);
        
        ColStart = ColStart+HColNum(L);
    end
  
    
    decoderData(qn0_1>=0) = 0;
    decoderData(qn0_1<0) = 1;
    if(mod(H*decoderData',2)==0)
        break;
    end
end


