%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdm_fading_sim4.m  %%%%%%%%%
%%%%%%%%%  data:2020年11月17日  author:飞蓬大将军 %%%%%%%%%%


%%%%%程序说明
%%%多径衰落信道下的OFDM传输
%%%调制方式：QPSK
%%%信道编码方式：无
%%%导频方式：梳状类型

%%%sim系列说明
%%%sim2:接收端不做捕获和同步，也不做信道估计，假设知道信道条件，后续版本考虑信道估计
%%%sim3:发送端插入梳状导频、接收端采用LS和MMSE信道估计算法
%%%sim4:尝试解决sim3中遗留问题，并将LS均衡方法、MMSE均衡方法与完美均衡(即接收端知道h，直接fft来做均衡)

%%%%    仿真环境
%软件版本：MATLAB R2019a

%********************** 程序主体 ************%

%%%%%%%%%%%%%%%%%%%%%  参数设置   %%%%%%%%%%%%%%%%%%%
para = 48;   %Number of parallel channel to transmit
fftlen = 64;  %FFT length
noc = 64;    %Number of carrier 
nd = 1;  %Number of information OFDM symbol for one loop
ml = 2;   %Modulation：QPSK
sr = 250000;  %Symbol rate 符号速率
br = sr.*ml;  %Bit rate per carrier
gilen = 16; %length of guard interval 

%%%%%%%%%%比特信噪比设置
%%%先设置大的，程序运行正确后，再设更多的信噪比设置，以画误码率曲线
% ebn0_temp = 80;
ebn0_temp = 1:1:30;
ber_fading_ls_linear = zeros(1,length(ebn0_temp));
ber_fading_ls_spline = zeros(1,length(ebn0_temp));
ber_fading_mmse = zeros(1,length(ebn0_temp));
ber_fading_per = zeros(1,length(ebn0_temp));

%%%%%%%导频信息
Nps = 4;  %导频间隔，起始
B_pilot = 1;  %导频起始位置，第1号子载波
Np = fftlen/Nps; %导频数量

for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 


%%%%%%%%%%%%%%%%%%%%%%%%%%  Fading initialization %%%%%%%%%%%%%%

PowerdB=[-1 -8 -17 -21 -25]; % 信道抽头功率特性
Delay=[0 3 5 6 8];          % 信道时延,示例
% Delay=[0 3 5 56 78];          % 信道时延
Power=10.^(PowerdB/10);     % 信道抽头功率特性 '线性'
Ntap=length(PowerdB);       % 信道抽头数
Lch=Delay(end)+1;           % 信道长度


%%%%%%%%%%%%%%%%%%%%%  主循环 %%%%%%%%%%%%%

nloop = 20000; %Number of sumulation loops 

noe2_ls_linear_temp = 0; %Number of error data  of LS_linear
noe2_ls_spline_temp = 0; %Number of error data of  LS_spline
noe2_mmse_temp = 0; %Number of error data of MMSE
noe2_per_temp = 0; %Number of error data of perfect compensation


nod = 0; %Number of transmitted data 

for iii = 1:nloop
    %%%%%%%%%%%%%%%%%  发射机  %%%%%%%%%%%%%%%%%%%
   seldata = rand(1,para*nd*ml)>0.5;  %串行数据

%     seldata = ones(1,para*nd*ml);  %用于调试程序
    paradata = reshape(seldata,para,nd*ml); %串并转换

    [ich,qch] = qpskmod(paradata,para,nd,ml); %调制
    kmod = 1/sqrt(2);
    ich1 = ich.*kmod;
    qch1 = qch.*kmod;
    ch1 = ich1 + qch1*1j;
    
    %%%%%%%%%%%将导频进行插入
%     Xp = 2*(randn(1,Np)>0)-1;    % Pilot sequence generation
%     Xp = 2*randi([0 1],1,Np) - 1;
    Xp = 2*ones(Np,nd) - 1;
    %%%方式一：
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

      %%%%%方式二：可以先把导频位置和数据位置算出来，比如kk1，kk2，相应放入
    
    %%%%%%%%%%%% IFFT %%%%%%%%%%

    y = ifft(X);
    ich2 = real(y);
    qch2 = imag(y);
    
     %%%观察信号功率
%     spow2 = sum(sum(ich2.^2+ qch2.^2))/nd./para;
    
    
    %%%%%%%% 添加保护间隔 %%%%%%%%
    [ich3,qch3] = giins(ich2,qch2,fftlen,gilen,nd);
    fftlen2 = fftlen + gilen;
    
    
    %****************  Attenuation Calculation **************
%     %%%方式一：
%     spow = sum(ich3.^2+ qch3.^2)/nd./para;
%     attn = 0.5*spow*sr/br*10.^(-ebn0/10);
%     attn = sqrt(attn);
   
    
    %%%方式二：
%     snr = ebn0 + 10*log10(2);
%     attn = sqrt(10.^(-snr/10)*spow/2);
    %以上两种方式表达是一样的
    
    %%%方式三：
    
    spow4 = sum(ich3.^2+ qch3.^2)/nd./fftlen2;
    esn0 = ebn0 + 10*log10(para/fftlen2) + 10*log10(ml);
    attn2 = 0.5*spow4*10.^(-esn0/10);
    attn2 = sqrt(attn2);
    
 
    
%     aaa = 1; %用于调试

    %***************  衰落信道 Fading channel ***************%
    channel = (randn(1,Ntap) + 1j * randn(1,Ntap)).*sqrt(Power/2);
    h = zeros(1,Lch);
    h(Delay+1) = channel;
    y = conv(ich3 + 1j*qch3,h);
    ifade = real(y(:,1:length(ich3)));
    qfade = imag(y(:,1:length(ich3)));
    spow5 = sum(ifade.^2+ qfade.^2)/nd./fftlen2;
    
%     %%%方式四：
%     esn0 = ebn0 + 10*log10(para/fftlen2) + 10*log10(ml);
%     attn2 = 0.5*spow5*10.^(-esn0/10);
%     attn2 = sqrt(attn2);
%     
    %%%%%不管哪种方式，以attn传入
    attn =  attn2;
    
    
    %***********************  接收机 *******************%
    %%%%%%%%%%假设已经完美同步上
    %%%%%%%%%% AWGN addition %%%%%%%%%
    [ich4,qch4] = comb(ifade,qfade,attn);

    %%%%%%%% 去掉保护间隔 %%%%%%%%
    [ich5,qch5] = girem(ich4,qch4,fftlen2,gilen,nd);
    
    
    %%%%%%%%%%%%%%  FFT %%%%%%%%
    rx = ich5 + qch5.*1i;
    ry = fft(rx);
%     ich6 = real(ry);
%     qch6 = imag(ry);

    %%%%%%%  信道估计  %%%%
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
    
   
    
    %%%%%%%%  信道均衡 %%%%
    %%%注意A的共轭转置和转置的区别，前者是A'，后者是A.'
    %%%%采用导频估计出的h，进行均衡
    ry_ls_linear_temp = ry./(H_est_ls_linear.');
    ry_ls_spline_temp = ry./(H_est_ls_spline.');
    ry_mmse_temp = ry./(H_est_mmse.');
    
    %%%假设接收端完美知道h，进行均衡
    H = fft([h,zeros(1,fftlen-Lch)].');
    ry_per_temp = ry./H; 
    
    
    
    
    %%%%%%%%%% 去除导频
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
    
    %%%%%%%%%%%%%% 并串转换 %%%%%%%%%
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

%************************* 画误码率曲线进行对比 ******************%
ebn0_temp = 1:1:30;
rayleign_one_path_theory = ber_temp(ebn0_temp); 
semilogy(ebn0_temp,rayleign_one_path_theory,'-*',ebn0_temp,ber_fading_ls_linear,'-^',ebn0_temp,ber_fading_ls_spline,'->',ebn0_temp,ber_fading_mmse,'-<',ebn0_temp,ber_fading_per,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('多径衰落信道下误码率仿真曲线');
legend('理论曲线','lslinear实验曲线','lsspline实验曲线','mmse实验曲线','完美均衡实验曲线');
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
%%%实验记录 2020年11月14日
%含已知导频的OFDM经过多径衰落信道的误码率仿真
%OFDM中引入保护间隔，是一种冗余信息，因此相比于理论误码率曲线有10*log（160/128）=0.969dB的损失
%OFDM信道估计引入了发送和接收端都已知的导频信号，这也是冗余，所以设计SNR与Eb/N0的换算
%本次实验还缺少相应的理论分析



