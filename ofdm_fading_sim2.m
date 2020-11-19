%%%%%%%%%%%%%%%%%%%%% OFDM���� %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdm_fading_sim2.m  %%%%%%%%%
%%%%%%%%%  data:2020��10��16��  author:����󽫾� %%%%%%%%%%


%%%%%����˵��
%%%�ྶ˥���ŵ��µ�OFDM����
%%%���Ʒ�ʽ��QPSK
%%%���뷽ʽ����
%%%���ն˲��������ͬ����Ҳ�����ŵ����ƣ�����֪���ŵ������������汾�����ŵ�����

%%%%    ���滷��
%����汾��MATLAB R2019a

%********************** �������� ************%

%%%%%%%%%%%%%%%%%%%%%  ��������   %%%%%%%%%%%%%%%%%%%
para = 128;   %Number of parallel channel to transmit
fftlen = 128;  %FFT length
noc = 128;    %Number of carrier 
nd = 6;  %Number of information OFDM symbol for one loop
ml = 2;   %Modulation��QPSK
sr = 250000;  %Symbol rate ��������
br = sr.*ml;  %Bit rate per carrier
gilen = 16; %length of guard interval 
ebn0_temp = 0:2:20;
ber_fading = zeros(1,length(ebn0_temp));
for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 


%%%%%%%%%%%%%%%%%%%%%%%%%%  Fading initialization %%%%%%%%%%%%%%

PowerdB=[0 -8 -17 -21 -25]; % �ŵ���ͷ��������
Delay=[0 3 5 6 8];          % �ŵ�ʱ��,ʾ��
% Delay=[0 3 5 56 78];          % �ŵ�ʱ��
Power=10.^(PowerdB/10);     % �ŵ���ͷ�������� '����'
Ntap=length(PowerdB);       % �ŵ���ͷ��
Lch=Delay(end)+1;           % �ŵ�����


%%%%%%%%%%%%%%%%%%%%%  ��ѭ�� %%%%%%%%%%%%%

nloop = 10000; %Number of sumulation loops 
noe = 0;   %Number of error data  
nod = 0;   %Number of transmitted data 
eop = 0;   %Number of error packet  
nop = 0;   %Number of transmitted packet


for iii = 1:nloop
    %%%%%%%%%%%%%%%%%  �����  %%%%%%%%%%%%%%%%%%%
   seldata = rand(1,para*nd*ml)>0.5;  %��������

%     seldata = ones(1,para*nd*ml); 
    paradata = reshape(seldata,para,nd*ml); %����ת��

    [ich,qch] = qpskmod(paradata,para,nd,ml); %����
    kmod = 1/sqrt(2);
    ich1 = ich.*kmod;
    qch1 = qch.*kmod;
    %%%%%%%%%%%% IFFT %%%%%%%%%%
    x = ich1 + qch1 *1j;
    
%     spow1 = sum(sum(ich1.^2+ qch1.^2))/nd./para;
    
    y = ifft(x);
    ich2 = real(y);
    qch2 = imag(y);
    
%     spow2 = sum(sum(ich2.^2+ qch2.^2))/nd./para;
    
    
    %%%%%%%% ��ӱ������ %%%%%%%%
    [ich3,qch3] = giins(ich2,qch2,fftlen,gilen,nd);
    fftlen2 = fftlen + gilen;
    
    
    %****************  Attenuation Calculation **************
    %%%��ʽһ��
    spow = sum(ich3.^2+ qch3.^2)/nd./para;
    attn = 0.5*spow*sr/br*10.^(-ebn0/10);
    attn = sqrt(attn);
    
    %%%��ʽ����
%     snr = ebn0 + 10*log10(2);
%     attn = sqrt(10.^(-snr/10)*spow/2);
    %�������ַ�ʽ�����һ����
    
    
    %***************  ˥���ŵ� Fading channel ***************%
    channel = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power/2);
    h = zeros(1,Lch);
    h(Delay+1) = channel;
    y = conv(ich3 + 1j*qch3,h);
    ifade = real(y(:,1:length(ich3)));
    qfade = imag(y(:,1:length(ich3)));
    
    %***********************  ���ջ� *******************%
    %%%%%%%%%% AWGN addition %%%%%%%%%
    [ich4,qch4] = comb(ifade,qfade,attn);

    %%%%%%%% ȥ��������� %%%%%%%%
    [ich5,qch5] = girem(ich4,qch4,fftlen2,gilen,nd);
    
    
    %%%%%%%%%%%%%%  FFT %%%%%%%%
    rx = ich5 + qch5.*1i;
    ry = fft(rx);
%     ich6 = real(ry);
%     qch6 = imag(ry);
    
    %%%%%%%%  �ŵ����� %%%%
    %%%ע��A�Ĺ���ת�ú�ת�õ�����ǰ����A'��������A.'
    H = fft([h,zeros(1,fftlen-Lch)].');
    
    for number = 1:nd
        ch6(:,number) = ry(:,number)./H; 
    end
    
    ich6 = real(ch6);
    qch6 = imag(ch6);
    
    %%%%%%%%%%%%%% demoluation %%%%%%%%%%%%%%
    ich7 = ich6./kmod;
    qch7 = qch6./kmod;
    demodata = qpskdemod(ich7,qch7,para,nd,ml);
    
    %%%%%%%%%%%%%% ����ת�� %%%%%%%%%
    demodata1 = reshape(demodata,1,para*nd*ml);
    
    %%%%%%%%%%%%%%% Bit Error Rate %%%%%%%%%%%
    noe2 = sum(abs(demodata1-seldata));
    nod2 = length(seldata);
    
    %%%%cumulative the number of error and data in noe and nod
    noe = noe + noe2;
    nod = nod + nod2;
     
     %%%����PER
    if noe2~=0
        eop = eop +1;

        
    end

    nop = nop + 1;
%     fprintf('%f\t%e\t%d\n',iii,noe2/nod2,eop);
    
%      fprintf('%f\t%e\t%d\n',iii,noe2/nod2,eop);
%     
    
end
%**************   output result *************
per = eop/nop;
ber = noe/nod;

ber_fading(1,kkk) = ber;
% fprintf('%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
% fid = fopen('BERofdm.dat','a');
% fprintf(fid,'%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
% fclose(fid);


end

%************************* �����������߽��жԱ� ******************%
rayleign_one_path_theory = ber_temp(ebn0_temp+1); 
semilogy(ebn0_temp,rayleign_one_path_theory,'-*',ebn0_temp,ber_fading,'-+');
xlabel('���������');
ylabel('������');
title('�ྶ˥���ŵ��������ʷ�������');
legend('��������','ʵ������');
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
%�����OFDM�����ྶ˥���ŵ��������ʷ���
%OFDM�����뱣���������һ��������Ϣ��������������������������10*log��160/128��=0.969dB����ʧ
%�ŵ������У��Ѿ����ŵ�������Ӱ�첹����
%2020��11��11��
