%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        *** Communication System ***                     %
%                            Single-Carrier FDMA(P6)                      %
%                                  Project                                %
%                               February, 2023                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%         By: SALFORE Alazar Gebrehiwot & TURTOI Rafael Andrei            %
%                                                                         %
%*************************************************************************%
%                                                                         %
% Description:This matlab code implements the computation of symbol error %
%              rate and plots the figure generated between and SNR vs SER %
%              for a single carrier frequency devision multiple access.   %
%                                                                         %
%*************************************************************************%                                                                         

clear all ; clc
close all;
%% =====================Input Parameters Definition======================%%

N = 128;                % The size of the transmitter IFFT and the receiver FFT
M = 16;                 % Symbols per block
Q = N/M;                % The bandwidth spreading factor for IFDMA
Q_D = Q-1;              % The bandwidth spreeding factor for a DFDMA
cp = 32 ;               % CP length.
SP.subband = 0;              
modulation_type='PSK';  
order = 2;              % BPSK 
equalizer_type='ZFEQ';
SNR = 0:1:30;           % Simulated SNR range is from 0 dB to 30 dB.
N_simu = 10^4;          % the number of simulation 


channel = [1 10^(-9.7/20) 10^(-22.8/20)];   % Channels based on 3GPP TS 25.104
channel = channel/sqrt(sum(channel.^2));    % Normalize the channel.


% Frequency expression of the channel response.
   H = fft(channel,N );
%% Transmitter
%=================== symbol error initialization =========================%
   Ifdma_SER=zeros(length(SNR),1);
   Dfdma_SER=zeros(length(SNR),1);
   Lfdma_SER=zeros(length(SNR),1);
    
for n = 1:length(SNR)
% Initialize the error count.
   Ifdma_error_number = 0;
   Dfdma_error_number = 0;
   Lfdma_error_number = 0;
for k = 1:N_simu
%=================== Generate random data block===========================%
% Input Generation
   input_symbol=randi([0, (order-1)] ,1,M);

%==============================Modulatation===============================%
% modulation_type=='PSK';
   input_signal = pskmod(input_symbol,order);

%========================Tranform to frequency domain=====================%
% Apply FFT operation
   input_signal_fft=fft(input_signal,M);

%===========================Subcarrier Mapping============================%
% Initialize subcarrier mappings 
   Ifdma_mapping = zeros(1,N);
   Dfdma_mapping = zeros(1,N);
   Lfdma_mapping = zeros(1,N);

% Applying Mapping
   Ifdma_mapping(1+SP.subband:Q:N) = input_signal_fft;
if Q ==1
   Dfdma_mapping=Ifdma_mapping;
else 
   Dfdma_mapping(1+SP.subband:Q_D:Q_D*M) = input_signal_fft; 

end 
   Lfdma_mapping([1:M]+M*SP.subband) = input_signal_fft;

%========================Tranform back to time domain=====================%
% Apply IFFT operation
   Ifdma_IFFT = ifft(Ifdma_mapping,N);
   Dfdma_IFFT = ifft(Dfdma_mapping,N);
   Lfdma_IFFT = ifft(Lfdma_mapping,N);

%==========================Add a cyclic prefix============================%
% Cyclic Prefix Addition
   Ifdma_cyclic = [Ifdma_IFFT(N-cp+1:N) Ifdma_IFFT ];
   Dfdma_cyclic = [Dfdma_IFFT(N-cp+1:N) Dfdma_IFFT ];
   Lfdma_cyclic = [Lfdma_IFFT(N-cp+1:N) Lfdma_IFFT ];
%% CHANNEL 

% multi-path channel
   Ifdma_signal=filter(channel, 1, Ifdma_cyclic);
   Dfdma_signal=filter(channel, 1, Dfdma_cyclic);
   Lfdma_signal=filter(channel, 1, Lfdma_cyclic);

% Generate AWGN 
   Noise  = (randn(1, N+cp)+1i*randn(1, N+cp))/sqrt(2); % N+cp by considering the added cyclic symbols 
   noisePower = 10^(-SNR(n)/10);

% Add AWGN to the transmitted signal.
   Ifdma_signal =Ifdma_signal+sqrt(noisePower/M)*Noise;
   Dfdma_signal =Dfdma_signal+sqrt(noisePower/M)*Noise;
   Lfdma_signal =Lfdma_signal +sqrt(noisePower/M)*Noise;
%% Receiver 
    
%=========================Remove the cyclic prefix========================%   
% Removing cyclic prefix
   Ifdma_received =  Ifdma_signal (cp+1:N+cp);
   Dfdma_received = Dfdma_signal(cp+1:N+cp);
   Lfdma_received = Lfdma_signal (cp+1:N+cp);

%========================Tranform to frequency domain=====================%  
% Applying FFT Operation
   Ifdma_received = fft(Ifdma_received,N);
   Dfdma_received = fft(Dfdma_received,N);
   Lfdma_received = fft(Lfdma_received,N);

%===========================Subcarrier De-Mapping=========================%
% Applying demaping
   Ifdma_received = Ifdma_received(1+SP.subband:Q:N);
if Q ==1
   Dfdma_received=Ifdma_received ;
else 
   Dfdma_received = Dfdma_received(1+SP.subband:Q_D:Q_D*M);   
end 
   
   Lfdma_received = Lfdma_received ([1:M]+M*SP.subband);
    
% channel response of the subcarriers
   H_Ifdma=H(1+SP.subband:Q:N);

if Q ==1
   H_Dfdma=  H_Ifdma ;
else 
   H_Dfdma=H(1+SP.subband:Q_D:Q_D*M);
end    
   H_Lfdma=H([1:M]+M*SP.subband);
   
%====================Perform frequency equalization=======================%
      
   Ifdma_received = Ifdma_received./H_Ifdma;
   Dfdma_received = Dfdma_received./H_Dfdma;
   Lfdma_received = Lfdma_received./H_Lfdma;

%========================Tranform back to time domain=====================%
% Applying  IFFT
   Ifdma_received = ifft(Ifdma_received);
   Dfdma_received = ifft(Dfdma_received);
   Lfdma_received = ifft(Lfdma_received);

%=============================De-Modulatation=============================%
% Symbol detection   
   Ifdma_symbol = pskdemod( Ifdma_received , order);
   Dfdma_symbol = pskdemod( Dfdma_received , order);
   Lfdma_symbol = pskdemod( Lfdma_received , order);
      
%% Error Calculation 
    
% Number of correctly received symbol 
   Ifdma_correct_symbol = find((input_symbol-Ifdma_symbol) == 0);
   Dfdma_correct_symbol = find((input_symbol-Dfdma_symbol) == 0);
   Lfdma_correct_symbol = find((input_symbol-Lfdma_symbol) == 0);
   
% The number of errors 
   Ifdma_error_number = Ifdma_error_number +(M-length(Ifdma_correct_symbol));
   Dfdma_error_number = Dfdma_error_number +(M-length(Dfdma_correct_symbol));
   Lfdma_error_number = Lfdma_error_number +(M-length(Lfdma_correct_symbol));
   
end 
% Calculating symbol error rate 
   Ifdma_SER(n,:) = Ifdma_error_number/ (M*N_simu);
   Dfdma_SER(n,:) = Dfdma_error_number/ (M*N_simu);
   Lfdma_SER(n,:) = Lfdma_error_number/ (M*N_simu);
end


%% plotting

% figure  
   semilogy(SNR,Ifdma_SER, '-*',SNR, Dfdma_SER, '-o',SNR, Lfdma_SER, '-d')

   first_zero1=find(Ifdma_SER==0,1,'first');

   Ifdma_SER=Ifdma_SER(1:first_zero1);

   first_zero2=find(Dfdma_SER==0,1,'first');
   Dfdma_SER=Dfdma_SER(1:first_zero2);

   first_zero3=find(Lfdma_SER==0,1,'first');
   Lfdma_SER=Lfdma_SER(1:first_zero3);

   Y_min=10^-8;

   semilogy(SNR(1:first_zero1),max(Y_min,Ifdma_SER),'LineWidth', 2)
   hold on, 

   semilogy(SNR(1:first_zero2), max(Y_min, Dfdma_SER),'LineWidth', 2)
   hold on, 

   semilogy(SNR(1:first_zero3), max(Y_min,Lfdma_SER),'LineWidth', 2)
   hold off

   legend(sprintf('IFDMA %K', Q),sprintf('DFDMA %K', Q),sprintf('LFDMA  %r', Q),'FontSize',10);

   title(['SER Diagram for SC-FDMA',' and input block size ',num2str(M)]);

   xlabel('Signal to Noise Ratio in [dB]', 'FontSize',10 );
   ylabel('Symbol Error rate', 'FontSize',10 );

   set(gca,'FontSize',12);
   axis( [0 30 Y_min 1]);



