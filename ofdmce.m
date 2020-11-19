%%%%%%%%%%%%%%%%%%%%% OFDM���� %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdmce.m  %%%%%%%%%
%%%%%%%%%  data:2020��11��7��  author:����󽫾� %%%%%%%%%%


%%%%%����˵��
%%%�е�Ƶ���Ƶ�OFDM����
%%% GI ChannelEsimation GI data GI data GI data GI data GI data GI data��1 channel estimation symbol 6 symbols��

%%%%    ���滷��
%����汾��MATLAB R2019a

%%%    simϵ��˵��

%********************** �������� ************%

%%%%%%%%%%%%%%%%%%%%%  ��������   %%%%%%%%%%%%%%%%%%%
para = 52;   %Number of parallel channel to transmit
fftlen = 64;  %FFT length
noc = 53;    %Number of carrier 
nd = 6;  %Number of information OFDM symbol for one loop
knd = 1; %Number of known channel estimation (CE) OFDM symbol�����ڵ�Ƶ���Ƶķ�����
ml = 2;   %Modulation��QPSK
sr = 250000;  %Symbol rate ��������
br = sr.*ml;  %Bit rate per carrier
gilen = 16; %length of guard interval 
%ebn0 = 3;  %Eb/No

ebn0_temp = 3:1:10;
for kkk = 1:length(ebn0_temp)
    
ebn0 = ebn0_temp(kkk);  %Eb/No 

%*******************  fading initialization **************%

%time resolution
tstp = 1/sr/(fftlen+gilen);

%ÿ�����ĵ���ʱ��
itau = [0];

%ÿ������ƽ������
dlvl = [0];
n0 = [6];
th1 = [0.0];


itnd0 = nd*(fftlen+gilen)*10;%Number of fading counter to skip
itnd1 = [1000]; %Initial value of fading counter
now1 = 1; %Number of direct wave + Number of delayed wave
fd = 150;%Maximum Doppler frequency

%the variable flat
flat = 0;
%%%%% flat = 1, ֻ�з���˥��
%%%%% flat = 0, �з��Ⱥ���λ˥��

%*******************  ������ **************%

nloop = 10000; %Number of sumulation loops 
noe = 0;   %Number of error data  
nod = 0;   %Number of transmitted data 
eop = 0;   %Number of error packet  
nop = 0;   %Number of transmitted packet


for iii = 1:nloop
    %%%%%%%%%%%%%%%%%  �����  %%%%%%%%%%%%%%%%%%%
    seldata = rand(1,para*nd*ml)>0.5;  %��������

    paradata = reshape(seldata,para,nd*ml); %����ת��

    [ich,qch] = qpskmod(paradata,para,nd,ml); %����
    kmod = 1/sqrt(2);
    ich = ich.*kmod;
    qch = qch.*kmod;
    
%   spow = sum(sum(ich1.^2+ qch1.^2))/nd./para;
  
%  %**************** ���뵼Ƶ ******** 
    %CE data generation 
%     kndata = zeros(1,fftlen);
%     kndata0 = 2.*(rand(1,52)<0.5) - 1;
%     kndata(1,2:27) = kndata0(1,1:26);
%     kndata(1,39:64) = kndata0(1,27:52);
%     ceich = kndata;
    ceich = ones(1,64);  %I·��Ƶ   
    ceqch = zeros(1,64); %Q·��Ƶ
    %**************  data mapping ***********%
    [ich1_temp,qch1_temp] = crmapping(ich,qch,fftlen,nd);
    
    ich1 = [ceich.' ich1_temp];
    qch1 = [ceqch.' qch1_temp];
    
    spow_flag1 = sum(sum(ich1.^2+ qch1.^2))/nd./para;
    %%%%%%%%%%%% IFFT %%%%%%%%%%
    x = ich1 + qch1 *1j;
   
    
    y = ifft(x);
    ich2 = real(y);
    qch2 = imag(y);
    
    spow_flag2 = sum(sum(ich2.^2+ qch2.^2))/nd./para;
    
    %��ifft֮������Ϊԭ���ģ�1/fftlen��������fft֮������Ϊԭ����fftlen��
    %��spow_flag1/spow_flag2 = 64
    
    %%%%%%%% ��ӱ������ %%%%%%%%
    [ich3,qch3] = giins(ich2,qch2,fftlen,gilen,nd+knd);
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%  ���ջ�  *************%%%%%%
    %%%%%%%%%%% AWGN addition *******************
    [ich4,qch4] = comb(ifade,qfade,attn);
    

 %**************  1 path rayleigh Fading perfect compensation  ***********%
    % If you would like to simulate performance under perfect
    % compensation,please remove '*'
%     ifade2 = 1./ramp.*(rcos(1,:).*ich4 + rsin(1,:).*qch4);
%     qfade2 = 1./ramp.*(-rsin(1,:).*ich4 + rcos(1,:).*qch4);
   
    %%%%%%%% ȥ��������� %%%%%%%%
%     [ich5,qch5] = girem(ifade2,qfade2,fftlen2,gilen,nd);
   [ich5,qch5] = girem(ich4,qch4,fftlen2,gilen,nd+knd);
    
    %%%%%%%%%%%%%%  FFT %%%%%%%%
    rx = ich5 + qch5.*1i;
    ry = fft(rx);
    ich6 = real(ry);
    qch6 = imag(ry);
    
    %*********************** Fading compensation by CE symbol **********%
   
    %��Ƶ�����Ƿ��Ͷˡ����ն˶���֪
    ice0 = ich1(:,knd); 
    qce0 = qch1(:,knd);
    
    ice1 = ich6(:,knd);
    qce1 = qch6(:,knd);
    
    
     %************** �Լ���д **********************%
    iv = (ice1.*ice0 + qce1.*qce0)./(ice0.^2 + qce0.^2);
    qv = (ice0.*qce1 - ice1.*qce0)./(ice0.^2 + qce0.^2);
    
    
   
%     iv1 = real(1./(ice0.^2 + qce0.^2).*(ice0 + 1i* qce0).*(ice1 - 1i* qce1));
%     qv1 = imag(1./(ice0.^2 + qce0.^2).*(ice0 + 1i* qce0).*(ice1 - 1i* qce1));
    %%%%�������ַ�ʽ��Ч
    ieqv1 = [iv iv iv iv iv iv iv]; 
    qeqv1 = [qv qv qv qv qv qv qv]; 
    
    [sizex,sizey] = size(ich6);
    icompen = zeros(sizex,sizey);
    qcompen = zeros(sizex,sizey);
    
%     icompen(:,1) = (iv.*ich6(:,1) + qv.*qch6(:,1))./(iv.^2 + qv.^2 );
%     qcompen(:,1) = (-qv.*ich6(:,1) + iv.*qch6(:,1))./(iv.^2 + qv.^2 );
%     aaa = 1;
    for t = 1:knd + nd
        icompen(:,t) = (iv.*ich6(:,t) + qv.*qch6(:,t))./(iv.^2 + qv.^2 );
        qcompen(:,t) = (-qv.*ich6(:,t) + iv.*qch6(:,t))./(iv.^2 + qv.^2 );
    end
    
%     %*************** �ο��� ***************%
%     
%     iv1 = real(1./(ice1.^2 + qce1.^2).*(ice0 + 1i* qce0).*(ice1 - 1i* qce1));
%     qv1 = imag(1./(ice1.^2 + qce1.^2).*(ice0 + 1i* qce0).*(ice1 - 1i* qce1));
%     
%     ieqv1 = [iv iv iv iv iv iv iv]; 
%     qeqv1 = [qv qv qv qv qv qv qv]; 
%      
%     icompen = real((ich6+1i*qch6).*( ieqv1 + 1i* qeqv1));
%     qcompen = imag((ich6+1i*qch6).*( ieqv1 + 1i* qeqv1));
%     
    ich7 = icompen(:,knd+1:knd+nd);
    qch7 = qcompen(:,knd+1:knd+nd);
    
    %**************  data mapping ***********%
    [ich8,qch8] = crdemapping(ich7,qch7,fftlen,nd);
    
   
    %%%%%%%%%%%%%% demoluation %%%%%%%%%%%%%%
    ich9 = ich8./kmod;
    qch9 = qch8./kmod;
    demodata = qpskdemod(ich9,qch9,para,nd,ml);
    
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
%***********   1������˥���µ�OFDM����BPSK���Ƶ�������ֵ **********%

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


%************************    ʵ��ֵ       ***************%
% % % % % % % %  EbN0(dB)      ������         �����        ѭ������
% % % % % % % % 3.000000	2.282389e-01	9.877000e-01	10000	150
% % % % % % % % 4.000000	2.061040e-01	9.723000e-01	10000	150
% % % % % % % % 5.000000	1.844702e-01	9.473000e-01	10000	150
% % % % % % % % 6.000000	1.634806e-01	9.138000e-01	10000	150
% % % % % % % % 7.000000	1.440835e-01	8.711000e-01	10000	150
% % % % % % % % 8.000000	1.254591e-01	8.200000e-01	10000	150
% % % % % % % % 9.000000	1.084558e-01	7.649000e-01	10000	150
% % % % % % % % 10.000000	9.311346e-02	7.113000e-01	10000	150

