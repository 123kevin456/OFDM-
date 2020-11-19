%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdmda.m  %%%%%%%%%
%%%%%%%%%  data:2020年10月19日  author:飞蓬大将军 %%%%%%%%%%


%%%%% 未完成 

%%%%%程序说明
%%% GI data GI data GI data GI data GI data GI data（6 symbols）
 
%%%%    仿真环境
%软件版本：MATLAB R2019a

%%%    



%********************** 程序主体 ************%

%%%%%%%%%%%%%%%%%%%%%  参数设置   %%%%%%%%%%%%%%%%%%%
para = 52;   %Number of parallel channel to transmit
fftlen = 64;  %FFT length
noc = 53;    %Number of carrier 
nd = 6;  %Number of information OFDM symbol for one loop
ml = 2;   %Modulation：QPSK
sr = 250000;  %Symbol rate 符号速率
br = sr.*ml;  %Bit rate per carrier
gilen = 16; %length of guard interval 
%ebn0 = 3;  %Eb/No

ebn0_temp = 3:1:10;
for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 

%*******************  fading initialization **************%

%time resolution
tstp = 1/sr/(fftlen+gilen);

%每个径的到达时间
itau = [0];

%每个径的平均功率
dlvl = [0];


n0 = [6];
th1 = [0.0];

%Number of fading counter to skip
itnd0 = nd*(fftlen+gilen)*10;

%Initial value of fading counter
itnd1 = [1000];
now1 = 1;

%Maximum Doppler frequency
fd = 150;

%the variable flat
flat = 0;
%%%%% flat = 1, 只有幅度衰落
%%%%% flat = 0, 有幅度和相位衰落



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
    ich = ich.*kmod;
    qch = qch.*kmod;
    
    %**************  data mapping ***********%
    [ich1,qch1] = crmapping(ich,qch,fftlen,nd);
    
    %%%%%%%%%%%% IFFT %%%%%%%%%%
    x = ich1 + qch1 *1j;
   
    
    y = ifft(x);
    ich2 = real(y);
    qch2 = imag(y);
    
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
    
 %***************   Fading channel ***************%
    [ifade,qfade,ramp,rcos,rsin] = sefade(ich3,qch3,itau,dlvl,th1,n0,itnd1,now1,length(ich3),tstp,fd,flat);
    itnd1 = itnd1 +itnd0;
%     ich3 = ifade;
%     qch3 = qfade;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%  接收机  *************%%%%%%
    %%%%%%%%%%% AWGN addition *******************
    [ich4,qch4] = comb(ifade,qfade,attn);
    

 %**************  1 path rayleigh Fading perfect compensation  ***********%
    % If you would like to simulate performance under perfect
    % compensation,please remove '*'
    ifade2 = 1./ramp.*(rcos(1,:).*ich4 + rsin(1,:).*qch4);
    qfade2 = 1./ramp.*(-rsin(1,:).*ich4 + rcos(1,:).*qch4);
  
    %%%%%%%% 去掉保护间隔 %%%%%%%%
    [ich5,qch5] = girem(ifade2,qfade2,fftlen2,gilen,nd);
 
    
    %%%%%%%%%%%%%%  FFT %%%%%%%%
    rx = ich5 + qch5.*1i;
    ry = fft(rx);
    ich6 = real(ry);
    qch6 = imag(ry);
    
    %**************  data mapping ***********%
    [ich7,qch7] = crdemapping(ich6,qch6,fftlen,nd);
    
    
    
    %%%%%%%%%%%%%% demoluation %%%%%%%%%%%%%%
    ich7 = ich7./kmod;
    qch7 = qch7./kmod;
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
% ber = noe/nod;
ber(kkk) = noe/nod;
% fprintf('%f\t%e\t%e\t%d\t%d\n',ebn0,ber,per,nloop,fd);
fprintf('%f\t%e\t%e\t%d\t%d\n',ebn0,ber(kkk),per,nloop,fd);
fid = fopen('BERofdma.dat','a');
% fprintf(fid,'%f\t%e\t%e\t%d\t%d\t\n',ebn0,ber,per,nloop,fd);
fprintf(fid,'%f\t%e\t%e\t%d\t%d\t\n',ebn0,ber(kkk),per,nloop,fd);
fclose(fid);
end
%%%%%%%%%%%%  simualtion result %%%%%%%%%%%

    
    
    
    
   
    
    
    
    
    
    
    
  




