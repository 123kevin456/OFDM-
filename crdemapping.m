%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% crdemapping.m  %%%%%%%%%
%%%%%%%%%  data:2020年10月19日  author:飞蓬大将军 %%%%%%%%%%

%Function to set data on subcarrier

function [iout,qout] = crdemapping(idata,qdata,fftlen,nd)

iout(1:26,:) = idata(2:27,:);
iout(27:52,:) = idata(39:64,:);

qout(1:26,:) = qdata(2:27,:);
qout(27:52,:) = qdata(39:64,:);


end
