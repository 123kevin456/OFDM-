%%%%%%%%%%%%%%%%%%%%%%%%%%   ������  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%     comb.m      %%%%%%%%%
%%%%%%%%%  data:2020��10��15��  author:����󽫾� %%%%%%%%%%


%********************** �������� ************%
function [iout,qout] = comb(idata,qdata,attn)

iout1 = randn(1,length(idata)).*attn;
qout1 = randn(1,length(qdata)).*attn;


iout = idata + iout1 ;
qout = qdata + qout1;


end