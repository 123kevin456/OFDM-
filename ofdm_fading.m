%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdm_fading.m  %%%%%%%%%
%%%%%%%%%  data:2020年10月16日  author:飞蓬大将军 %%%%%%%%%%


%%%%%程序说明
%%%一径衰落信道下的OFDM传输

%%%%    仿真环境
%软件版本：MATLAB R2019a

%%%    sim系列说明


%********************** 程序主体 ************%

%%%%%%%%%%%%%%%%%%%%%  参数设置   %%%%%%%%%%%%%%%%%%%
para = 128;   %Number of parallel channel to transmit
fftlen = 128;  %FFT length
noc = 128;    %Number of carrier 
nd = 6;  %Number of information OFDM symbol for one loop
ml = 2;   %Modulation：QPSK
sr = 250000;  %Symbol rate 符号速率
br = sr.*ml;  %Bit rate per carrier
gilen = 32; %length of guard interval 

ebn0_temp = 3:1:10;
for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 


%%%%%%%%%%%%%%%%%%%%%%%%%%  Fading initialization %%%%%%%%%%%%%%


%time resolution
tstp = 1/sr/(fftlen+gilen);

%每个径的到达时间
itau = [0];  %Delay time for each multipath fading

%每个径的平均功率
dlvl = [0]; %Attenuation level for each multipath fading


n0 = [6];  %Number of waves in order to generate each multipath fading
th1 = [0.0]; %Initialized phase for each multipath fading

%Number of fading counter to skip
itnd0 = nd*(fftlen+gilen)*10;

%Initial value of fading counter
itnd1 = [1000];
now1 = 1;

%Maximum Doppler frequency
fd = 320;

%the variable flat
flat = 1;
%%%%% flat = 1, 只有幅度衰落
%%%%% flat = 0, 有幅度和相位衰落

%%%%%%%%%%%%%%%%%%%%%  主循环 %%%%%%%%%%%%%

nloop = 10000; %Number of sumulation loops 
noe = 0;   %Number of error data  
nod = 0;   %Number of transmitted data 
eop = 0;   %Number of error packet  
nop = 0;   %Number of transmitted packet


for iii = 1:nloop
    %%%%%%%%%%%%%%%%%  发射机  %%%%%%%%%%%%%%%%%%%
    seldata = rand(1,para*nd*ml)>0.5;  %串行数据

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
    spow = sum(ich3.^2+ qch3.^2)/nd./para;
    attn = 0.5*spow*sr/br*10.^(-ebn0/10);
    attn = sqrt(attn);

    %***************  衰落信道 Fading channel ***************%
    [ifade,qfade] = sefade(ich3,qch3,itau,dlvl,th1,n0,itnd1,now1,length(ich3),tstp,fd,flat);
    
    %%%%
    itnd1 = itnd1 +itnd0;
    
    %***********************  接收机 *******************%
    %%%%%%%%%% AWGN addition %%%%%%%%%
    [ich4,qch4] = comb(ifade,qfade,attn);
     

   
      %%%%%%%% 去掉保护间隔 %%%%%%%%
    [ich5,qch5] = girem(ich4,qch4,fftlen2,gilen,nd);
    
    
    
    %%%%%%%%%%%%%%  FFT %%%%%%%%
    rx = ich5 + qch5.*1i;
    ry = fft(rx);
    ich6 = real(ry);
    qch6 = imag(ry);
    
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
%     else
%         eop = eop;
        
    end

    nop = nop + 1;
%     fprintf('%f\t%e\t%d\n',iii,noe2/nod2,eop);
    
%      fprintf('%f\t%e\t%d\n',iii,noe2/nod2,eop);
%     
    
end
%**************   output result *************
per = eop/nop;
ber = noe/nod;

fprintf('%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fid = fopen('BERofdm.dat','a');
fprintf(fid,'%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fclose(fid);


end

%************************* 画误码率曲线进行对比 ******************%
ebn0 = 3:1:10;
rayleign_one_path_theory = [0.0919131757263162,0.0771369160563911,0.0641826854495229,0.0529988839256387,0.0434744067460626,0.0354590676278380,0.0287823671004334,0.0232687053772039]; 
rayleigh_flat_0 = [1.379502e-01,1.215418e-01,1.059700e-01,9.183802e-02,7.866901e-02,6.675632e-02,5.598379e-02,4.649043e-02]; 
semilogy(ebn0,rayleign_one_path_theory,'-*',ebn0,rayleigh_flat_0,'-+');
xlabel('比特信噪比');
ylabel('误码率');
title('不同信噪比下误码率仿真曲线');
legend('理论曲线','实验曲线');
grid on;




%***********   1径瑞利衰落下的OFDM采用BPSK调制的误码率值 **********%

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

%%%%%%%%%%%%%%%%%%%%%%%      实验值          **************%

% flat = 0

%%%%%%%%%%%%%     EbN0(dB)      误码率         误包率        循环次数
% % % % %      3.000000	     4.999334e-01	1.000000e+00	10000	
% % % % %      4.000000	     5.001447e-01	1.000000e+00	10000	
% % % % %      5.000000	     4.998994e-01	1.000000e+00	10000	
% % % % %      6.000000	     5.001079e-01	1.000000e+00	10000	
% % % % %      7.000000	     5.001012e-01	1.000000e+00	10000	
% % % % %      8.000000	     5.000331e-01	1.000000e+00	10000	
% % % % %      9.000000	     5.001562e-01	1.000000e+00	10000	
% % % % %     10.000000	     4.999725e-01	1.000000e+00	10000	

% flat = 1
%%%%%%%%%%%%%     EbN0(dB)      误码率         误包率        循环次数
% % % % % %      3.000000	1.379502e-01	9.422000e-01	10000	
% % % % % %      4.000000	1.215418e-01	9.060000e-01	10000	
% % % % % %      5.000000	1.059700e-01	8.622000e-01	10000	
% % % % % %      6.000000	9.183802e-02	8.075000e-01	10000	
% % % % % %      7.000000	7.866901e-02	7.536000e-01	10000	
% % % % % %     8.000000	6.675632e-02	6.931000e-01	10000	
% % % % % %     9.000000	5.598379e-02	6.359000e-01	10000	
% % % % % %    10.000000	4.649043e-02	5.725000e-01	10000	


















