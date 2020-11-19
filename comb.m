%%%%%%%%%%%%%%%%%%%%%%%%%%   加噪声  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%     comb.m      %%%%%%%%%
%%%%%%%%%  data:2020年10月15日  author:飞蓬大将军 %%%%%%%%%%


%********************** 程序主体 ************%
function [iout,qout] = comb(idata,qdata,attn)

iout1 = randn(1,length(idata)).*attn;
qout1 = randn(1,length(qdata)).*attn;


iout = idata + iout1 ;
qout = qdata + qout1;


end