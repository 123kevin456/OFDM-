%%%%%%%%%%%%%%%%%%%%% OFDM���� %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% crmapping.m  %%%%%%%%%
%%%%%%%%%  data:2020��10��19��  author:����󽫾� %%%%%%%%%%

%Function to set data on subcarrier

function [iout,qout] = crmapping(idata,qdata,fftlen,nd)
iout = zeros(fftlen,nd);
qout = zeros(fftlen,nd);

iout(2:27,:) = idata(1:26,:);
iout(39:64,:) = idata(27:52,:);

qout(2:27,:) = qdata(1:26,:);
qout(39:64,:) = qdata(27:52,:);


end

