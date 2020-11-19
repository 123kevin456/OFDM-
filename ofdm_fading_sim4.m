%%%%%%%%%%%%%%%%%%%%% OFDM���� %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdm_fading_sim4.m  %%%%%%%%%
%%%%%%%%%  data:2020��11��17��  author:����󽫾� %%%%%%%%%%


%%%%%����˵��
%%%�ྶ˥���ŵ��µ�OFDM����
%%%���Ʒ�ʽ��QPSK
%%%�ŵ����뷽ʽ����
%%%��Ƶ��ʽ����״����

%%%simϵ��˵��
%%%sim2:���ն˲��������ͬ����Ҳ�����ŵ����ƣ�����֪���ŵ������������汾�����ŵ�����
%%%sim3:���Ͷ˲�����״��Ƶ�����ն˲���LS��MMSE�ŵ������㷨
%%%sim4:���Խ��sim3���������⣬����LS���ⷽ����MMSE���ⷽ������������(�����ն�֪��h��ֱ��fft��������)

%%%%    ���滷��
%����汾��MATLAB R2019a

%********************** �������� ************%

%%%%%%%%%%%%%%%%%%%%%  ��������   %%%%%%%%%%%%%%%%%%%
para = 48;   %Number of parallel channel to transmit
fftlen = 64;  %FFT length
noc = 64;    %Number of carrier 
nd = 1;  %Number of information OFDM symbol for one loop
ml = 2;   %Modulation��QPSK
sr = 250000;  %Symbol rate ��������
br = sr.*ml;  %Bit rate per carrier
gilen = 16; %length of guard interval 

%%%%%%%%%%�������������
%%%�����ô�ģ�����������ȷ������������������ã��Ի�����������
% ebn0_temp = 80;
ebn0_temp = 1:1:30;
ber_fading_ls_linear = zeros(1,length(ebn0_temp));
ber_fading_ls_spline = zeros(1,length(ebn0_temp));
ber_fading_mmse = zeros(1,length(ebn0_temp));
ber_fading_per = zeros(1,length(ebn0_temp));

%%%%%%%��Ƶ��Ϣ
Nps = 4;  %��Ƶ�������ʼ
B_pilot = 1;  %��Ƶ��ʼλ�ã���1�����ز�
Np = fftlen/Nps; %��Ƶ����

for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 


%%%%%%%%%%%%%%%%%%%%%%%%%%  Fading initialization %%%%%%%%%%%%%%

PowerdB=[-1 -8 -17 -21 -25]; % �ŵ���ͷ��������
Delay=[0 3 5 6 8];          % �ŵ�ʱ��,ʾ��
% Delay=[0 3 5 56 78];          % �ŵ�ʱ��
Power=10.^(PowerdB/10);     % �ŵ���ͷ�������� '����'
Ntap=length(PowerdB);       % �ŵ���ͷ��
Lch=Delay(end)+1;           % �ŵ�����


%%%%%%%%%%%%%%%%%%%%%  ��ѭ�� %%%%%%%%%%%%%

nloop = 20000; %Number of sumulation loops 

noe2_ls_linear_temp = 0; %Number of error data  of LS_linear
noe2_ls_spline_temp = 0; %Number of error data of  LS_spline
noe2_mmse_temp = 0; %Number of error data of MMSE
noe2_per_temp = 0; %Number of error data of perfect compensation


nod = 0; %Number of transmitted data 

for iii = 1:nloop
    %%%%%%%%%%%%%%%%%  �����  %%%%%%%%%%%%%%%%%%%
   seldata = rand(1,para*nd*ml)>0.5;  %��������

%     seldata = ones(1,para*nd*ml);  %���ڵ��Գ���
    paradata = reshape(seldata,para,nd*ml); %����ת��

    [ich,qch] = qpskmod(paradata,para,nd,ml); %����
    kmod = 1/sqrt(2);
    ich1 = ich.*kmod;
    qch1 = qch.*kmod;
    ch1 = ich1 + qch1*1j;
    
    %%%%%%%%%%%����Ƶ���в���
%     Xp = 2*(randn(1,Np)>0)-1;    % Pilot sequence generation
%     Xp = 2*randi([0 1],1,Np) - 1;
    Xp = 2*ones(Np,nd) - 1;
    %%%��ʽһ��
    X = zeros(fftlen,nd);
    ip = 0;    
    pilot_loc = [];
    for k=1:fftlen
       if mod(k,Nps)==1
             temp = floor(k/Nps)+1;
             X(k,:) = Xp(temp,:);
%            X(k) = Xp(floor(k/Nps)+1); 
             pilot_loc = [pilot_loc k]; 
             ip = ip+1;
         else
             X(k,:) = ch1(k-ip,:);
       end
    end

      %%%%%��ʽ���������Ȱѵ�Ƶλ�ú�����λ�������������kk1��kk2����Ӧ����
    
    %%%%%%%%%%%% IFFT %%%%%%%%%%

    y = ifft(X);
    ich2 = real(y);
    qch2 = imag(y);
    
     %%%�۲��źŹ���
%     spow2 = sum(sum(ich2.^2+ qch2.^2))/nd./para;
    
    
    %%%%%%%% ��ӱ������ %%%%%%%%
    [ich3,qch3] = giins(ich2,qch2,fftlen,gilen,nd);
    fftlen2 = fftlen + gilen;
    
    
    %****************  Attenuation Calculation **************
%     %%%��ʽһ��
%     spow = sum(ich3.^2+ qch3.^2)/nd./para;
%     attn = 0.5*spow*sr/br*10.^(-ebn0/10);
%     attn = sqrt(attn);
   
    
    %%%��ʽ����
%     snr = ebn0 + 10*log10(2);
%     attn = sqrt(10.^(-snr/10)*spow/2);
    %�������ַ�ʽ�����һ����
    
    %%%��ʽ����
    
    spow4 = sum(ich3.^2+ qch3.^2)/nd./fftlen2;
    esn0 = ebn0 + 10*log10(para/fftlen2) + 10*log10(ml);
    attn2 = 0.5*spow4*10.^(-esn0/10);
    attn2 = sqrt(attn2);
    
 
    
%     aaa = 1; %���ڵ���

    %***************  ˥���ŵ� Fading channel ***************%
    channel = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power/2);
    h = zeros(1,Lch);
    h(Delay+1) = channel;
    y = conv(ich3 + 1j*qch3,h);
    ifade = real(y(:,1:length(ich3)));
    qfade = imag(y(:,1:length(ich3)));
    spow5 = sum(ifade.^2+ qfade.^2)/nd./fftlen2;
    
%     %%%��ʽ�ģ�
%     esn0 = ebn0 + 10*log10(para/fftlen2) + 10*log10(ml);
%     attn2 = 0.5*spow5*10.^(-esn0/10);
%     attn2 = sqrt(attn2);
%     
    %%%%%�������ַ�ʽ����attn����
    attn =  attn2;
    
    
    %***********************  ���ջ� *******************%
    %%%%%%%%%%�����Ѿ�����ͬ����
    %%%%%%%%%% AWGN addition %%%%%%%%%
    [ich4,qch4] = comb(ifade,qfade,attn);

    %%%%%%%% ȥ��������� %%%%%%%%
    [ich5,qch5] = girem(ich4,qch4,fftlen2,gilen,nd);
    
    
    %%%%%%%%%%%%%%  FFT %%%%%%%%
    rx = ich5 + qch5.*1i;
    ry = fft(rx);
%     ich6 = real(ry);
%     qch6 = imag(ry);

    %%%%%%%  �ŵ�����  %%%%
    for m=1:3
         if m==1 
             H_est_ls_linear = LS_CE(ry.',Xp.',pilot_loc,fftlen,Nps,'linear');
%              method='LS-linear'; % LS estimation with linear interpolation
         elseif m==2
             H_est_ls_spline = LS_CE(ry.',Xp.',pilot_loc,fftlen,Nps,'spline'); 
%              method='LS-spline'; % LS estimation with spline interpolation
         else
%              H_est_mmse = MMSE_CE(ry.',Xp.',pilot_loc,fftlen,Nps,h,SNR); 
              H_est_mmse = MMSE_CE(ry.',Xp.',pilot_loc,fftlen,Nps,h,ebn0); 
%              method='MMSE'; % MMSE estimation
         end
    end % end for count'
    
   
    
    %%%%%%%%  �ŵ����� %%%%
    %%%ע��A�Ĺ���ת�ú�ת�õ�����ǰ����A'��������A.'
    %%%%���õ�Ƶ���Ƴ���h�����о���
    ry_ls_linear_temp = ry./(H_est_ls_linear.');
    ry_ls_spline_temp = ry./(H_est_ls_spline.');
    ry_mmse_temp = ry./(H_est_mmse.');
    
    %%%������ն�����֪��h�����о���
    H = fft([h,zeros(1,fftlen-Lch)].');
    ry_per_temp = ry./H; 
    
    
    
    
    %%%%%%%%%% ȥ����Ƶ
    ip = 0;
    for k=1:fftlen
       if mod(k,Nps)==1
%              temp = floor(k/Nps)+1;
%              X(k,:) = Xp(temp,:);
% %            X(k) = Xp(floor(k/Nps)+1); 
%              pilot_loc = [pilot_loc k]; 
             ip = ip+1;
       else
             ry_ls_linear(k-ip,:) = ry_ls_linear_temp(k,:);
             ry_ls_spline(k-ip,:) = ry_ls_spline_temp(k,:);
             ry_mmse(k-ip,:) = ry_mmse_temp(k,:);
             ry_per(k-ip,:) = ry_per_temp(k,:);
       end
    end
    
    
    %%%%%%%%%%%%%% demoluation %%%%%%%%%%%%%
    demodata_ls_linear = qpskdemod(real(ry_ls_linear)./kmod,imag(ry_ls_linear)./kmod,para,nd,ml);
    demodata_ls_spline= qpskdemod(real(ry_ls_spline)./kmod,imag(ry_ls_spline)./kmod,para,nd,ml);
    demodata_mmse = qpskdemod(real(ry_mmse)./kmod,imag(ry_mmse)./kmod,para,nd,ml);
    demodata_per = qpskdemod(real(ry_per)./kmod,imag(ry_per)./kmod,para,nd,ml);
    
    %%%%%%%%%%%%%% ����ת�� %%%%%%%%%
    demodata1_ls_linear = reshape(demodata_ls_linear,1,para*nd*ml);
    demodata1_ls_spline = reshape(demodata_ls_spline,1,para*nd*ml);
    demodata1_mmse  = reshape(demodata_mmse,1,para*nd*ml);
    demodata1_per = reshape(demodata_per,1,para*nd*ml);
    
    
    %%%%%%%%%%%%%%% Bit Error Rate %%%%%%%%%%%
    noe2_ls_linear = sum(abs(demodata1_ls_linear-seldata));
    noe2_ls_spline = sum(abs(demodata1_ls_spline-seldata));
    noe2_mmse = sum(abs(demodata1_mmse-seldata));
    noe2_per = sum(abs(demodata1_per-seldata));
    nod2 = length(seldata);
    
    %%%%cumulative the number of error and data in noe and nod
    noe2_ls_linear_temp = noe2_ls_linear_temp + noe2_ls_linear; %Number of error data  of LS_linear
    noe2_ls_spline_temp = noe2_ls_spline_temp + noe2_ls_spline ; %Number of error data of  LS_spline
    noe2_mmse_temp = noe2_mmse_temp + noe2_mmse; %Number of error data of MMSE
    noe2_per_temp = noe2_per_temp + noe2_per;
    nod = nod + nod2;
        
    
end
%**************   output result *************
ber_ls_linear_temp = noe2_ls_linear_temp/nod;
ber_ls_spline_temp = noe2_ls_spline_temp/nod;
ber_mmse_temp = noe2_mmse_temp/nod;
ber_per_temp = noe2_per_temp/nod;


ber_fading_ls_linear(1,kkk) = ber_ls_linear_temp ;
ber_fading_ls_spline(1,kkk) = ber_ls_spline_temp;
ber_fading_mmse(1,kkk) = ber_mmse_temp;
ber_fading_per(1,kkk) = ber_per_temp;


end

%************************* �����������߽��жԱ� ******************%
ebn0_temp = 1:1:30;
rayleign_one_path_theory = ber_temp(ebn0_temp); 
semilogy(ebn0_temp,rayleign_one_path_theory,'-*',ebn0_temp,ber_fading_ls_linear,'-^',ebn0_temp,ber_fading_ls_spline,'->',ebn0_temp,ber_fading_mmse,'-<',ebn0_temp,ber_fading_per,'-+');
xlabel('���������');
ylabel('������');
title('�ྶ˥���ŵ��������ʷ�������');
legend('��������','lslinearʵ������','lssplineʵ������','mmseʵ������','��������ʵ������');
grid on;
%***********   �ྶ����˥���µ�OFDM����BPSK���Ƶ�������ֵ **********%

%%%%%%%%%%%%%%%%%%%%%%%      ����ֵ          **************%
%%%%%%%%%%%%%     EbN0(dB)      ������        
%%%%%%%%%%%%%       3        0.125000000000000
%%%%%%%%%%%%%       4        0.100000000000000
%%%%%%%%%%%%%       5        0.0833333333333333
%%%%%%%%%%%%%       6        0.0714285714285715
%%%%%%%%%%%%%       7        0.0625000000000000
%%%%%%%%%%%%%       8        0.0555555555555556
%%%%%%%%%%%%%       9        0.0500000000000000
%%%%%%%%%%%%%      10        0.0454545454545455

%%%%%%%%%%%%%%%%%   ����    %%%%%%%%%%%%%%%%%%%%
%%%ʵ���¼ 2020��11��14��
%����֪��Ƶ��OFDM�����ྶ˥���ŵ��������ʷ���
%OFDM�����뱣���������һ��������Ϣ��������������������������10*log��160/128��=0.969dB����ʧ
%OFDM�ŵ����������˷��ͺͽ��ն˶���֪�ĵ�Ƶ�źţ���Ҳ�����࣬�������SNR��Eb/N0�Ļ���
%����ʵ�黹ȱ����Ӧ�����۷���



