function [ iter,decoderData ] = ldpcdecoderminsum(H,HRowNum,HColNum,receiveSignal,MAX_ITER_NUM,NORM_FACTOR)
%LDPCC decode algorithm , based on the Min-Sum algorithm
% H: check matrix
% HRowNum,HColNum: index generate from check matrix
% receiveSignal: received soft information from channel
% MAT_INTER_NUM: maximum iterations
% HRowNum,HColNum generating by the following codes
% 	[r_mark,c_mark] = find(H~=0);
% 	HColNum = sum(H);
% 	HRowNum = cell(1,size(H,1));
% 	for rowH = 1:size(H,1)
%		 HRowNum{rowH} = find(r_mark==rowH);
% 	end

[~,N] = size(H);    

vl = receiveSignal;
decoderData = zeros(1,N);

uml = zeros(1,sum(HColNum));
vml = uml;
ColStart = 1;
for L=1:length(HColNum)
    vml(ColStart:ColStart+HColNum(L)-1) = vl(L);
    ColStart = ColStart+HColNum(L);
end

for iter=1:MAX_ITER_NUM
    %check nodes information process
    for L=1:length(HRowNum)
        L_col = HRowNum{L};
        vmltemp = vml(L_col);
        vmlMark = ones(size(vmltemp));
        vmlMark(vmltemp<0) = -1;
        vmlMark = prod(vmlMark);
        minvml = sort(abs(vmltemp));
        for L_col_i = 1:length(L_col)
            if minvml(1)==abs(vmltemp(L_col_i))
                if vmltemp(L_col_i)<0
                    vmltemp(L_col_i) = -vmlMark*minvml(2);
                else
                    vmltemp(L_col_i) = vmlMark*minvml(2);
                end
            else
                if vmltemp(L_col_i)<0
                    vmltemp(L_col_i) = -vmlMark*minvml(1);
                else
                    vmltemp(L_col_i) = vmlMark*minvml(1);
                end
            end
        end
        uml(L_col) = NORM_FACTOR*vmltemp;
    end
    %variable nodes information process
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
    % decision decoding
    decoderData(qn0_1>=0) = 0;
    decoderData(qn0_1<0) = 1;
    if(mod(H*decoderData',2)==0)
        break;
    end
end


