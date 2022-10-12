clear all
close all
%% SYMULATION FLAG
Rayleigh = false; %Set to true to apply Rayleigh channel, false is AWGN
waterfilling = false; %Set to true to apply waterfilling on Rayleigh channel
numOfBlocks = 30; %Number of simulated blocks: it takes time!!
%% Parameters definition
global W;
W= 2*10^6;
global Tc;
Tc=1/W;
global Nc;
Nc = 130;
global Lmax;
Lmax=12;
global h;
%% Channel creation and Frequency Response
h = getChannelRealization();
bin = W/(3*Nc);
for ii=1:4
    h(ii)=normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
end
h(10)=normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
h(Lmax)=normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));

%Normalization
h = channelNormalization(h);
f = 0:bin:(W+2*bin);
H = zeros(1,length(f));
for ii=1:length(f)
    for jj=1:Lmax
        H(ii)=H(ii)+h(jj)*exp(-1i*2*pi*f(ii)/W*jj);
    end
end
ch_fft  = fft(h);
%% Data initialization: BPSK without noise with simbols +1 and -1
d=(2*randi([0,1],1,Nc))-1; %Creation of symbols in frequency
x = rayleighSig(d);
y=channelApply(x);
y=y(Lmax:end);
res = (1/sqrt(Nc))*fft(y);
%% STEP 3-Channel inversion and error estimation
res = res./ch_fft;
errore = sum((abs(d-res)).^2);
%% STEP 4-QPSK modulation
P = 10;
SNRdB = 2:2:14;
SNR = 10.^(SNRdB/10);
%Symbols definition, with unitary average power
symbolsQPSK = [(1+1i),(1-1i),(-1+1i),(-1+-1i)]; %modulation type=2
symbolsQPSK = symbolsQPSK./sqrt(mean(abs(symbolsQPSK).^2));
bitTableQPSK = [0,1,2,3];

symbolsBPSK = [-1,1]; %modulation type=1
symbolsBPSK = symbolsBPSK./sqrt(mean(abs(symbolsBPSK).^2));
bitTableBPSK = [0,1];

symbols16QAM = [3*(-1+1i),-1+3*1i,1+3*1i,3*(1+1i),-3+1i,-1+1i,1+1i,3+1i,-3-1i,-1-1i,1-1i,3-1i,-3-3*1i,-1-3*1i,1-3*1i,3-3*1i];
symbols16QAM = symbols16QAM./sqrt(mean(abs(symbols16QAM).^2));
bitTable16QAM= [0,4,12,8,1,5,13,9,3,7,15,11,2,6,14,10];%modulation type=4

messageSymbolsLength = numOfBlocks*Nc;
if(Rayleigh)
   rep_channel_fft = zeros(numOfBlocks,Nc);
   rep_channel_fft(1,:)=ch_fft;
   for ii=2:(numOfBlocks)
        y = fft(channelNormalization(getChannelRealization()));       
        rep_channel_fft(ii,:) = y;
   end
end
%Pbit definition
Pbit = zeros(1,length(SNR));
confInt= zeros(2,length(SNR));
for ii=1:length(SNR)
    modulationType = ones(1,messageSymbolsLength);
    N0= P/SNR(ii);
    newNoiseVector = normrnd(0,sqrt((N0/2)),[1,messageSymbolsLength]) +1i*normrnd(0,sqrt((N0/2)),[1,messageSymbolsLength]); 
    if(~waterfilling)
        messageIndexes=randi([1 length(symbolsQPSK)],1,messageSymbolsLength);
        symbolsSent = sqrt(P)*symbolsQPSK(messageIndexes);
        bitsSent=reshape(decimalToBinaryVector(bitTableQPSK(messageIndexes),2).',1,[]);
        powerUsed = P*ones(1,messageSymbolsLength);
        modulationType = 2*modulationType;
    end
    if(Rayleigh)
        rep_channel_fft1 = reshape(rep_channel_fft.',1,[]);
        newNoiseVector = newNoiseVector./rep_channel_fft1;
    end
    %encoding symbols waterfilling
    if(Rayleigh && waterfilling)
        powerUsed = ones(1,messageSymbolsLength);
        symbolsSent = ones(1,messageSymbolsLength);
        messageIndex = ones(1,messageSymbolsLength);
        %power definition
        tmp = 1;
            for jj=1:numOfBlocks
                [pow,lda]=applyWaterfilling(rep_channel_fft(jj,:),P,N0);
                powerUsed(tmp:tmp+Nc-1) = pow;
                tmp = tmp+Nc;
            end 
         %symbols creation
         bitsSent = zeros(1,messageSymbolsLength*4);
         nxtBit = 1;
         aborted_symbols=0;
            for jj=1:messageSymbolsLength
               if(powerUsed(jj)==0)
                   modulationType(jj) = 0;
                   aborted_symbols = aborted_symbols+1;
                   symbolsSent(jj)=0;
               elseif(powerUsed(jj)<P/2)
                       messageIndex = randi([1 length(symbolsBPSK)]);
                       numberSent = bitTableBPSK(messageIndex);
                       symbolsSent(jj) = sqrt(powerUsed(jj))*symbolsBPSK(messageIndex);
                       modulationType(jj) = 1;                                     
               elseif (powerUsed(jj)<P+0.1*P || SNR(ii)<14)
                        messageIndex = randi([1 length(symbolsQPSK)]);
                        numberSent = bitTableQPSK(messageIndex);
                        symbolsSent(jj) = sqrt(powerUsed(jj))*symbolsQPSK(messageIndex);
                        modulationType(jj) = 2;
                               else
                                messageIndex = randi([1 length(symbols16QAM)]);
                                numberSent = bitTable16QAM(messageIndex);
                                symbolsSent(jj) = sqrt(powerUsed(jj))*symbols16QAM(messageIndex);
                                modulationType(jj) = 4;  
               end                      
               if(modulationType(jj) ~=0)
                    bitsSent(nxtBit:nxtBit+modulationType(jj)-1) = decimalToBinaryVector(numberSent,modulationType(jj));
               end
                nxtBit=nxtBit+modulationType(jj);
            end
    end
    %fine waterfilling
    zeroElemIndex = find(~modulationType);
    symbolsReceived=symbolsSent + newNoiseVector;
    symbolsReceived(zeroElemIndex)=0;
    symbolsSent(symbolsSent~=0);
    modulationType(modulationType~=0);
    totBit = sum(modulationType) ;
    
    if(waterfilling)
        bitsSent = bitsSent(1:totBit);
    end
    
    bitRx = ones(1,totBit);
    %decoding
    posBit = 1;
    for jj=1:messageSymbolsLength
        switch modulationType(jj)
            case 0
               nbit =0; 
            case 2
                symbols = symbolsQPSK;
                nbit =2;
                bitTable = bitTableQPSK;
            case 1
                symbols = symbolsBPSK;
                nbit =1;
                bitTable = bitTableBPSK;           
            case 4
                symbols = symbols16QAM;
                nbit =4;
                bitTable = bitTable16QAM;
            otherwise
                disp('Error')
        end
        if(modulationType(jj) ~=0)
            symbolsReceivedPos = decode(symbolsReceived(jj),symbols,powerUsed(jj));
            numberRx = bitTable(symbolsReceivedPos); %numero di rx bits
            bitRx(posBit:posBit+nbit-1) = decimalToBinaryVector(numberRx,nbit);
        end
        posBit = posBit+nbit;
    end
      
     %error rx
     errorReception = sum(abs(bitsSent-bitRx));
     nerrs = errorReception;    % Number of bit errors in simulation
     ntrials = length(bitsSent); % Number of trials in simulation
     level = 0.95;   % Confidence level
     [Pbit(ii),confInt(:,ii)] = berconfint(nerrs,ntrials,level);
    
end



PbitTeorico = qfunc(sqrt(SNR));
%Waterfilling
N0= P/SNR(5);
[powerprofile,lda] = applyWaterfilling(ch_fft,P,N0);

% PLOTTING SECTION
%% Rayleigh and Waterfilling-Rayleigh Pbit
%NOTE: to obtain the comparisons graphs, I saved the data of these graphs
%and then plotted on the same one. All results can be obtained, but
%singularly (not on the same graph).

%NOTE:this plot shows up only if the flag Rayleigh is set to true;
if(Rayleigh)
figure
hold on
grid on
plot(SNRdB,log(Pbit),'bo-');
plot(SNRdB,log(confInt(1,:)),"bx");
plot(SNRdB,log(confInt(2,:)),"bx");
title('BER with confidence interval as function of SNR');
xlabel("SNR[dB]");
ylabel("log(Pbit)");
legend('Pbit Rayleigh with confidence interval')
end
%% Gaussian Pbit
%NOTE:use this plot only if the flag Rayleigh is set to false;
if(~Rayleigh)
figure
hold on
grid on
plot(SNRdB(1:5),log(Pbit(1:5)),'ro-');
plot(SNRdB(1:5),log(PbitTeorico(1:5)),'bo-');
plot(SNRdB(1:5),log(confInt(1,(1:5))),"kx");
plot(SNRdB(1:5),log(confInt(2,(1:5))),"kx");
title('BER with confidence intervals over AWGN channel');
xlabel("SNR[dB]");
ylabel("log(Pbit)");
legend('Estimated QPSK BER with confidence intervals', 'Theoretical QPSK BER ')
end
%% Frequency response 
figure
hold on
plot(f,abs(H));
title('Frequency Response of the Baseband Channel');
xlabel("Hz");
xlim([0 W+2*bin])
f1 = (W/Nc)*(0:1:Nc-1);
plot(f1,abs(ch_fft),"o");
legend('H(f)','H(f=nW/Nc)')
%% Impulse response 
figure
grid on
l = 0:1:Lmax-1;
stem(l,abs(h(1:Lmax)),'o');
title('Normalized 6-taps channel');
xlabel("ℓ");
ylabel("|h_{ℓ}|");
%% Waterfilling
f1 = (W/Nc)*(0:1:Nc-1);
l = ones(length(f1))*1/lda;
figure
plot(f1,N0./abs(ch_fft).^2,"-o");
hold on
plot(f1,powerprofile,"-o");
plot(f1,l,"green");
legend('N_0/|H(f=nW/N_c)|^2','P*(f=nW/N_c)','1/\lambda');
xlim([0 W+2*bin])
ylim([0 1/lda+1]);
title('Optimal power allocation');
xlabel("Hz");
function decodedSymbolPosition = decode(rxSymbol,symbols,power)
    dist = zeros(1,length(symbols)); 
        for ii=1:length(symbols)
            dist(ii)=abs(rxSymbol-sqrt(power)*symbols(ii));
        end
        [minumum,index] = min(dist);
        decodedSymbolPosition=index;
end
function y=channelApply(signal)
global Lmax;
global h;
    y= zeros(1,length(signal));
        for ii=Lmax:length(signal)
            for jj=1:Lmax
                y(ii)=y(ii)+h(jj)*signal(ii-jj+1);
            end
        end
end
function sig = rayleighSig(freqSymbols)
global Lmax;
global Nc;
    d_idft = sqrt(Nc)*ifft(freqSymbols); %idft of symbols
    d_pre = d_idft(length(d_idft)-Lmax+2:length(d_idft)); %prefix creation
    
    sig = [d_pre,d_idft]; %Signal to be sent
end
function ch = getChannelRealization()
    global Lmax;
    global Nc;
    ch = zeros(1,Nc);
    for ii=1:4
        ch(ii)=normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
    end
    ch(10)=normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
    ch(Lmax)=normrnd(0,sqrt(0.5)) + 1i*normrnd(0,sqrt(0.5));
end
function c_norm = channelNormalization(h)
    c_norm = h/sqrt((sum(abs(h).^2)));
end
function [powerProfile,lda] = applyWaterfilling(h_freq,P,N0)
global Nc;
   h_sq=abs(h_freq).^2;
   ldaMin=0.001;
   fun=@(lda)abs((1/Nc)*sum(max((1/lda*ones(1,length(h_freq)) - N0./h_sq),0))-P);
   lda = fminsearch(fun,ldaMin);
   powerProfile = max((1/lda*ones(1,length(h_freq)) - N0./h_sq),0);
end