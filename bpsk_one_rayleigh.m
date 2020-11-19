


ebn0_temp = 1:1:40;
for kkk = 1:length(ebn0_temp)
    
% ebn0 = 10^(ebn0_temp(kkk)/10);  %Eb/No，换成十进制
ebn0 = (10^(ebn0_temp(kkk)/10))*(48/80);  %Eb/No，换成十进制
% ebn0 = ebn0_temp(kkk);  %Eb/No，dB形式
ebn0_temp2 = 1./(ebn0+1);
ebn0_temp3 = (1 -  ebn0_temp2)^(0.5);
ber_temp(1,kkk) = 0.5*(1-ebn0_temp3);

end
