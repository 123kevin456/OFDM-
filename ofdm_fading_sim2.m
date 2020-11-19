%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdm_fading_sim2.m  %%%%%%%%%
%%%%%%%%%  data:2020年10月16日  author:飞蓬大将军 %%%%%%%%%%


%%%%%程序说明
%%%多径衰落信道下的OFDM传输
%%%调制方式：QPSK
%%%编码方式：无
%%%接收端不做捕获和同步，也不做信道估计，假设知道信道条件，后续版本考虑信道估计

%%%%    仿真环境
%软件版本：MATLAB R2019a

%********************** 程序主体 ************%

%%%%%%%%%%%%%%%%%%%%%  参数设置   %%%%%%%%%%%%%%%%%%%
para = 128;   %Number of parallel channel to transmit
fftlen = 128;  %FFT length
noc = 128;    %Number of carrier 
nd = 6;  %Number of information OFDM symbol for one loop
ml = 2;   %Modulation：QPSK
sr = 250000;  %Symbol rate 符号速率
br = sr.*ml;  %Bit rate per carrier
gilen = 16; %length of guard interval 
ebn0_temp = 0:2:20;
ber_fading = zeros(1,length(ebn0_temp));
for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 


%%%%%%%%%%%%%%%%%%%%%%%%%%  Fading initialization %%%%%%%%%%%%%%

PowerdB=[0 -8 -17 -21 -25]; % 信道抽头功率特性
Delay=[0 3 5 6 8];          % 信道时延,示例
% Delay=[0 3 5 56 78];          % 信道时延
Power=10.^(PowerdB/10);     % 信道抽头功率特性 '线性'
Ntap=length(PowerdB);       % 信道抽头数
Lch=Delay(end)+1;           % 信道长度


%%%%%%%%%%%%%%%%%%%%%  主循环 %%%%%%%%%%%%%

nloop = 10000; %Number of sumulation loops 
noe = 0;   %Number of error data  
nod = 0;   %Number of transmitted data 
eop = 0;   %Number of error packet  
nop = 0;   %Number of transmitted packet


for iii = 1:nloop
    %%%%%%%%%%%%%%%%%  发射机  %%%%%%%%%%%%%%%%%%%
   seldata = rand(1,para*nd*ml)>0.5;  %串行数据

%     seldata = ones(1,para*nd*ml); 
    paradata = reshape(seldata,para,nd*ml); %串并转换

    [ich,qch] = qpskmod(paradata,para,nd,ml); %调制
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
    
    
    %%%%%%%% 添加保护间隔 %%%%%%%%
    [ich3,qch3] = giins(ich2,qch2,fftlen,gilen,nd);
    fftlen2 = fftlen + gilen;
    
    
    %****************  Attenuation Calculation **************
    %%%方式一：
    spow = sum(ich3.^2+ qch3.^2)/nd./para;
    attn = 0.5*spow*sr/br*10.^(-ebn0/10);
    attn = sqrt(attn);
    
    %%%方式二：
%     snr = ebn0 + 10*log10(2);
%     attn = sqrt(10.^(-snr/10)*spow/2);
    %以上两种方式表达是一样的
    
    
    %***************  衰落信道 Fading channel ***************%
    channel = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power/2);
    h = zeros(1,Lch);
    h(Delay+1) = channel;
    y = conv(ich3 + 1j*qch3,h);
    ifade = real(y(:,1:length(ich3)));
    qfade = imag(y(:,1:length(ich3)));
    
    %***********************  接收机 *******************%
    %%%%%%%%%% AWGN addition %%%%%%%%%
    [ich4,qch4] = comb(ifade,qfade,attn);

    %%%%%%%% 去掉保护间隔 %%%%%%%%
    [ich5,qch5] = girem(ich4,qch4,fftlen2,gilen,nd);
    
    
    %%%%%%%%%%%%%%  FFT %%%%%%%%
    rx = ich5 + qch5.*1i;
    ry = fft(rx);
%     ich6 = real(ry);
%     qch6 = imag(ry);
    
    %%%%%%%%  信道均衡 %%%%
    %%%注意A的共轭转置和转置的区别，前者是A'，后者是A.'
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
    
    %%%%%%%%%%%%%% 并串转换 %%%%%%%%%
    demodata1 = reshape(demodata,1,para*nd*ml);
    
    %%%%%%%%%%%%%%% Bit Error Rate %%%%%%%%%%%
    noe2 = sum(abs(demodata1-seldata));
    nod2 = length(seldata);
    
    %%%%cumulative the number of error and data in noe and nod
    noe = noe + noe2;
    nod = nod + nod2;
     
     %%%计算PER
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

%************************* 画误码率曲线进行对比 ******************%
rayleign_one_path_theory = ber_temp(ebn0_temp+1); 
semilogy(ebn0_temp,rayleign_one_path_theory,'-*',ebn0_temp,ber_fading,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('多径衰落信道下误码率仿真曲线');
legend('理论曲线','实验曲线');
grid on;
%***********   多径瑞利衰落下的OFDM采用BPSK调制的误码率值 **********%

%%%%%%%%%%%%%%%%%%%%%%%      理论值          **************%
%%%%%%%%%%%%%     EbN0(dB)      误码率        
%%%%%%%%%%%%%       3        0.125000000000000
%%%%%%%%%%%%%       4        0.100000000000000
%%%%%%%%%%%%%       5        0.0833333333333333
%%%%%%%%%%%%%       6        0.0714285714285715
%%%%%%%%%%%%%       7        0.0625000000000000
%%%%%%%%%%%%%       8        0.0555555555555556
%%%%%%%%%%%%%       9        0.0500000000000000
%%%%%%%%%%%%%      10        0.0454545454545455

%%%%%%%%%%%%%%%%%   结论    %%%%%%%%%%%%%%%%%%%%
%完成了OFDM经过多径衰落信道的误码率仿真
%OFDM中引入保护间隔，是一种冗余信息，因此相比于理论误码率曲线有10*log（160/128）=0.969dB的损失
%信道均衡中，已经将信道带来的影响补偿上
%2020年11月11日
