%% LDPCC encode and decode simulation
%
%    Autor: Cao Wenhui
%    Last Modify:2015-03-07
%    Runtime:MATLAB(R) 2014a
%

%% 
%% 
% clear everything 
clc;clear all;close all;

%%
% Simulation parameter setting
EbN0_dB = 1.5:0.2:2;                       
FRAMES_NUM = 10;                         
MAX_ITER_NUM = 250;                      
MAX_ERROR_FRAME = 200;                   
bitError = zeros(1,length(EbN0_dB));     
BER = bitError;                          
frameError = bitError;                  
iterNumTotal = zeros(1,length(EbN0_dB)); 
INFO_LENGTH = 1024;                      
RATE = 1/2;                              
SIZE_M = 512;                            

%%
% open(create) a file to save important data.
FILE_NAME = ['LDPC_CCSDS_' datestr(now,'yyyymmdd') '.txt'];
fid = fopen(FILE_NAME,'at+');
fprintf(fid,'date %s\n',datestr(now,'yyyymmdd'));
fprintf(fid,'%s\n',' LDPCC in CCSDS standard section3(2006)');
fprintf(fid,'Information Length = %d, ',INFO_LENGTH);
fprintf(fid,'Information Rate = %d, ',RATE);

%% 
% H represents check matrix;G represents generate matrix
H = ccsdscheckmatrix(SIZE_M,RATE);
G = ccsdsgeneratematrix(H,SIZE_M,RATE);
%% 
% preprocessing
[r_mark,c_mark] = find(H~=0);
HColNum = sum(H);
HRowNum = cell(1,size(H,1));
for rowH = 1:size(H,1)
    HRowNum{rowH} = find(r_mark==rowH);
end
%% BER Monte Carlo Simulation

%%
% BER in different EbN0
for nEbN0 = 1:length(EbN0_dB)
%%
% Monte Carlo Simulation
    for nF=1:FRAMES_NUM
        %%
        % encode
        message = randi([0 1],1,INFO_LENGTH);
        encodeData = mod(message*G,2);
        %%
        % modulate
        transmitSignal = 1 - 2*encodeData;   % 0-1;1--1
        transmitSignalPower = sqrt(var(transmitSignal));
        transmitSignal = transmitSignal/transmitSignalPower; %Normalization
        %%
        % AWGN Channel, the relationship between SNR and EbN0
        %%
        % 
        % $$SNR = \frac{{{E_b}}}{{{N_0}}} + 2{\log _{10}}^M + 2{\log _{10}}^{RATE}$$
        % 
        
        SNR_dB = EbN0_dB((nEbN0)) + 10*log10(2)+10*log10(RATE);
        SNR = 10^(SNR_dB/10);
        noise = randn(1,length(transmitSignal));
        noise = noise/sqrt(SNR);     
        receiveSignal = transmitSignal + noise;
        %%
        % punching
        receiveSignal(end-1*SIZE_M+1:end-0*SIZE_M) = 0;
        %%
        % decode
        [iterNum,recoverData] = ...      
            ldpcdecoderllr(H,HRowNum,HColNum,receiveSignal,SNR,MAX_ITER_NUM);
            %ldpcdecoderbp(H,HRowNum,HColNum,receiveSignal,SNR,MAX_ITER_NUM);
            %ldpcdecoderminsum(H,HRowNum,HColNum,receiveSignal,SNR,MAX_ITER_NUM);
        % output
        if(nEbN0==1 && nF==1)
            fprintf(fid,'decoding function is ldpcdecoderllr\n');
            fprintf(fid,'MAX-iteration = %d, ',MAX_ITER_NUM);
            fprintf(fid,'\n-------------------------------\n');
            fprintf(fid,' table name in Chinese\n');
            
            fprintf(fid,'EbN0\t T_Frame\t E_Frames\t E_Bits\t A_IterNums\t BER\t FER\n');
        end
        %% 
        % BER,FER and iterations
        bitError(nEbN0) = bitError(nEbN0) + ...
            sum(abs(message - recoverData(1:length(message))));
        frameError(nEbN0) = frameError(nEbN0) + ...
            (sum(abs(encodeData - recoverData(1:length(encodeData))))~=0);
        iterNumTotal(nEbN0) = iterNumTotal(nEbN0) + iterNum;
        %% 
        % 
         if ( frameError(nEbN0)>=MAX_ERROR_FRAME || nF==FRAMES_NUM)
            BER(nEbN0) = bitError(nEbN0)/nF/length(message);
            fprintf(fid,'%3.2g\t %5d\t %4d\t ',EbN0_dB(nEbN0),nF,frameError(nEbN0));
            fprintf(fid,'%8d\t %8.6g\t ',bitError(nEbN0),iterNumTotal(nEbN0)/nF);
            fprintf(fid,'%e\t %e\n',BER(nEbN0),frameError(nEbN0)/nF);       
             break;
         end
         if (mod(nF,100)==0)
            BER(nEbN0) = bitError(nEbN0)/nF/320;
            fprintf('\n-------------------------------\n');
            fprintf('Eb/No = %e \n',EbN0_dB(nEbN0));
            fprintf('Total Frames = %d, ',nF);
            fprintf('Error Frames = %d, ',frameError(nEbN0));
            fprintf('Bit error = %d, ',bitError(nEbN0));
            fprintf('MAX-iteration = %d, ',MAX_ITER_NUM);
            fprintf('Average-interation = %e, ',iterNumTotal(nEbN0)/nF);
            fprintf('BER = %g, ',BER(nEbN0));
            fprintf('FER = %g \n',frameError(nEbN0)/nF);            
         end
    end
end

%% simulation result
%%
fclose(fid);
readtable(FILE_NAME,'HeaderLines',6,'Delimiter','\t')
semilogy(EbN0_dB,BER,'r');
xlabel('Eb/N0(dB)'); ylabel('BER');
grid on;