%%%%%%%%%%%%%%%%%%%%%%%%%%   ȥ���������  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%     girem.m      %%%%%%%%%
%%%%%%%%%  data:2020��10��15��  author:����󽫾� %%%%%%%%%%

%********************** �������� ************%
function [iout,qout] = girem(idata,qdata,fftlen2,gilen,nd)

idata2 = reshape(idata,fftlen2,nd);
qdata2 = reshape(qdata,fftlen2,nd);

iout = idata2(gilen+1:fftlen2,:);
qout = qdata2(gilen+1:fftlen2,:);


end