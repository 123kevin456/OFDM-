%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdm.m  %%%%%%%%%
%%%%%%%%%  date:2020年10月11日  author:飞蓬大将军 %%%%%%%%%%


%%%%%程序说明
%%%高斯白噪声信道下OFDM传输

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
gilen = 32; %length of guard interval 


ebn0_temp = 3:1:10;
for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 
% ebn0 = 8;  %Eb/No

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
    
    spow1 = sum(sum(ich1.^2+ qch1.^2))/nd./para;
    
    y = ifft(x);
    ich2 = real(y);
    qch2 = imag(y);
    
    spow2 = sum(sum(ich2.^2+ qch2.^2))/nd./para;
    
    
    %%%%%%%% 添加保护间隔 %%%%%%%%
    [ich3,qch3] = giins(ich2,qch2,fftlen,gilen,nd);
    fftlen2 = fftlen + gilen;
    
%     figure(1);
%     plot(ich3,'-+');
%     hold on;
    
    
    
    %****************  Attenuation Calculation **************
    spow = sum(ich3.^2+ qch3.^2)/nd./para;
    attn = 0.5*spow*sr/br*10.^(-ebn0/10);
    attn = sqrt(attn);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%  接收机  *************%%%%%%
    %%%%%%%%%%% AWGN addition *******************
    [ich4,qch4] = comb(ich3,qch3,attn);
    

%     plot(ich4,'-*');
  
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
% ebn0 = 3:1:10;
% awgn_theory = [0.0228784075610853,0.0125008180407376,0.00595386714777866,0.00238829078093281,0.000772674815378444,0.000190907774075993,3.36272284196176e-05,3.87210821552205e-06];
% awgn_no_compensation = [3.698496e-02,2.254329e-02,1.226654e-02,5.823633e-03,2.305339e-03,7.492187e-04,1.757812e-04,3.170573e-05];
% rayleign_one_path_theory = [0.125000000000000,0.100000000000000,0.0833333333333333,0.0714285714285715,0.0625000000000000,0.0555555555555556,0.0500000000000000,0.0454545454545455]; 
% 
% semilogy(ebn0,awgn_theory,'-*',ebn0,awgn_no_compensation,'-+');
% xlabel('比特信噪比');
% ylabel('误码率');
% title('不同信噪比下误码率仿真曲线');
% legend('理论曲线','实验曲线');
% grid on;

%***********   高斯白噪声下的OFDM采用QPSK调制的误码率值 **********%

%%%%%%%%%%%%%%%%%%%%%%%      理论值          **************%
%%%%%%%%%%%%%     EbN0(dB)      误码率        
%%%%%%%%%%%%%       3        0.0228784075610853
%%%%%%%%%%%%%       4        0.0125008180407376
%%%%%%%%%%%%%       5        0.00595386714777866
%%%%%%%%%%%%%       6        0.00238829078093281
%%%%%%%%%%%%%       7        0.000772674815378444
%%%%%%%%%%%%%       8        0.000190907774075993
%%%%%%%%%%%%%       9        3.36272284196176e-05
%%%%%%%%%%%%%      10        3.87210821552205e-06

%%%%%%%%%%%%%%%%%%%%%        实验值        *******%%%%%%%%%%%%
%%%%%%%%%%%%%      EbN0(dB)      误码率         误包率        循环次数
% % % % % % % % % 3.000000	3.698496e-02	1.000000e+00	10000	
% % % % % % % % % 4.000000	2.254329e-02	1.000000e+00	10000	
% % % % % % % % % 5.000000	1.226654e-02	1.000000e+00	10000	
% % % % % % % % % 6.000000	5.823633e-03	9.999000e-01	10000	
% % % % % % % % % 7.000000	2.305339e-03	9.709000e-01	10000	
% % % % % % % % % 8.000000	7.492187e-04	6.837000e-01	10000	
% % % % % % % % % 9.000000	1.757812e-04	2.384000e-01	10000	
% % % % % % % % % 10.000000	3.170573e-05	4.780000e-02	10000	    
   
%%%%%%%%%%%%%%%%%   结论    %%%%%%%%%%%%%%%%%%%%
%完成了OFDM经过高斯白噪声信道的误码率仿真
%由于OFDM中引入保护间隔，是一种冗余信息，因此相比于理论误码率曲线有10*log（160/128）=0.969dB的损失
%对应到误码率曲线上，即本实验误码率曲线相比理论误码率曲线向右平移0.969dB
%2020年11月7日   
   
    
    
    
    
    
    
    
  




