%%%%%%%%%%%%%%%%%%%%%%%%%%   �ӱ������  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%     giins.m      %%%%%%%%%
%%%%%%%%%  data:2020��10��15��  author:����󽫾� %%%%%%%%%%

%********************** �������� ************%
function [iout,qout] = giins(idata,qdata,fftlen,gilen,nd)
idata1 = reshape(idata,fftlen,nd);
qdata1 = reshape(qdata,fftlen,nd);


idata2 = [idata1(fftlen - gilen +1:fftlen,:); idata1];
qdata2 = [qdata1(fftlen - gilen +1:fftlen,:); qdata1];

iout = reshape(idata2,1,(fftlen+gilen)*nd);
qout = reshape(qdata2,1,(fftlen+gilen)*nd);


end