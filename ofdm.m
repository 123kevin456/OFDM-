%%%%%%%%%%%%%%%%%%%%% OFDM���� %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% ofdm.m  %%%%%%%%%
%%%%%%%%%  date:2020��10��11��  author:����󽫾� %%%%%%%%%%


%%%%%����˵��
%%%��˹�������ŵ���OFDM����

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
    %%%%%%%%%%%%%%%%%  �����  %%%%%%%%%%%%%%%%%%%
    seldata = rand(1,para*nd*ml)>0.5;  %��������

    paradata = reshape(seldata,para,nd*ml); %����ת��

    [ich,qch] = qpskmod(paradata,para,nd,ml); %����
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
    
    
    %%%%%%%% ��ӱ������ %%%%%%%%
    [ich3,qch3] = giins(ich2,qch2,fftlen,gilen,nd);
    fftlen2 = fftlen + gilen;
    
%     figure(1);
%     plot(ich3,'-+');
%     hold on;
    
    
    
    %****************  Attenuation Calculation **************
    spow = sum(ich3.^2+ qch3.^2)/nd./para;
    attn = 0.5*spow*sr/br*10.^(-ebn0/10);
    attn = sqrt(attn);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%  ���ջ�  *************%%%%%%
    %%%%%%%%%%% AWGN addition *******************
    [ich4,qch4] = comb(ich3,qch3,attn);
    

%     plot(ich4,'-*');
  
    %%%%%%%% ȥ��������� %%%%%%%%
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
ber = noe/nod;

fprintf('%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fid = fopen('BERofdm.dat','a');
fprintf(fid,'%f\t%e\t%e\t%d\t\n',ebn0,ber,per,nloop);
fclose(fid);

end

%************************* �����������߽��жԱ� ******************%
% ebn0 = 3:1:10;
% awgn_theory = [0.0228784075610853,0.0125008180407376,0.00595386714777866,0.00238829078093281,0.000772674815378444,0.000190907774075993,3.36272284196176e-05,3.87210821552205e-06];
% awgn_no_compensation = [3.698496e-02,2.254329e-02,1.226654e-02,5.823633e-03,2.305339e-03,7.492187e-04,1.757812e-04,3.170573e-05];
% rayleign_one_path_theory = [0.125000000000000,0.100000000000000,0.0833333333333333,0.0714285714285715,0.0625000000000000,0.0555555555555556,0.0500000000000000,0.0454545454545455]; 
% 
% semilogy(ebn0,awgn_theory,'-*',ebn0,awgn_no_compensation,'-+');
% xlabel('���������');
% ylabel('������');
% title('��ͬ������������ʷ�������');
% legend('��������','ʵ������');
% grid on;

%***********   ��˹�������µ�OFDM����QPSK���Ƶ�������ֵ **********%

%%%%%%%%%%%%%%%%%%%%%%%      ����ֵ          **************%
%%%%%%%%%%%%%     EbN0(dB)      ������        
%%%%%%%%%%%%%       3        0.0228784075610853
%%%%%%%%%%%%%       4        0.0125008180407376
%%%%%%%%%%%%%       5        0.00595386714777866
%%%%%%%%%%%%%       6        0.00238829078093281
%%%%%%%%%%%%%       7        0.000772674815378444
%%%%%%%%%%%%%       8        0.000190907774075993
%%%%%%%%%%%%%       9        3.36272284196176e-05
%%%%%%%%%%%%%      10        3.87210821552205e-06

%%%%%%%%%%%%%%%%%%%%%        ʵ��ֵ        *******%%%%%%%%%%%%
%%%%%%%%%%%%%      EbN0(dB)      ������         �����        ѭ������
% % % % % % % % % 3.000000	3.698496e-02	1.000000e+00	10000	
% % % % % % % % % 4.000000	2.254329e-02	1.000000e+00	10000	
% % % % % % % % % 5.000000	1.226654e-02	1.000000e+00	10000	
% % % % % % % % % 6.000000	5.823633e-03	9.999000e-01	10000	
% % % % % % % % % 7.000000	2.305339e-03	9.709000e-01	10000	
% % % % % % % % % 8.000000	7.492187e-04	6.837000e-01	10000	
% % % % % % % % % 9.000000	1.757812e-04	2.384000e-01	10000	
% % % % % % % % % 10.000000	3.170573e-05	4.780000e-02	10000	    
   
%%%%%%%%%%%%%%%%%   ����    %%%%%%%%%%%%%%%%%%%%
%�����OFDM������˹�������ŵ��������ʷ���
%����OFDM�����뱣���������һ��������Ϣ��������������������������10*log��160/128��=0.969dB����ʧ
%��Ӧ�������������ϣ�����ʵ���������������������������������ƽ��0.969dB
%2020��11��7��   
   
    
    
    
    
    
    
    
  




