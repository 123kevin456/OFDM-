%%%%%%%%%%%%%%%%%%%%%%%%%%   去除保护间隔  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%     girem.m      %%%%%%%%%
%%%%%%%%%  data:2020年10月15日  author:飞蓬大将军 %%%%%%%%%%

%********************** 程序主体 ************%
function [iout,qout] = girem(idata,qdata,fftlen2,gilen,nd)

idata2 = reshape(idata,fftlen2,nd);
qdata2 = reshape(qdata,fftlen2,nd);

iout = idata2(gilen+1:fftlen2,:);
qout = qdata2(gilen+1:fftlen2,:);


end