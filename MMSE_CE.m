function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR)
%function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,ts,SNR)
% MMSE channel estimation function
% Inputs:
%       Y         = Frequency-domain received signal
%       Xp        = Pilot signal
%       pilot_loc = Pilot location
%       Nfft      = FFT size
%       Nps       = Pilot spacing
%       h         = Channel impulse response
%       ts        = Sampling time
%       SNR       = Signal-to-Noise Ratio[dB]
% output:
%      H_MMSE     = MMSE channel estimate

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

%H = fft(h,N);
snr = 10^(SNR*0.1);
Np=Nfft/Nps; 
k=1:Np; 
H_tilde = Y(1,pilot_loc(k))./Xp(k);  % LS estimate
k=0:length(h)-1; %k_ts = k*ts; 
hh = h*h'; 
tmp = h.*conj(h).*k; %tmp = h.*conj(h).*k_ts;
r = sum(tmp)/hh;    
r2 = tmp*k.'/hh; %r2 = tmp*k_ts.'/hh;
tau_rms = sqrt(r2-r^2);     % rms delay
df = 1/Nfft;  %1/(ts*Nfft);
j2pi_tau_df = j*2*pi*tau_rms*df;
K1 = repmat([0:Nfft-1].',1,Np); 
K2 = repmat([0:Np-1],Nfft,1);
rf = 1./(1+j2pi_tau_df*(K1-K2*Nps));
K3 = repmat([0:Np-1].',1,Np); 
K4 = repmat([0:Np-1],Np,1);
rf2 = 1./(1+j2pi_tau_df*Nps*(K3-K4));
Rhp = rf;
Rpp = rf2 + eye(length(H_tilde),length(H_tilde))/snr;
H_MMSE = transpose(Rhp*inv(Rpp)*H_tilde.');  % MMSE channel estimate