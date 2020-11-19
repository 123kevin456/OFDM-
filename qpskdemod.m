%%%%%%%%%%%%%%%%%%%%%%%%%%   qpsk解调  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   qpskdemod.m      %%%%%%%%%
%%%%%%%%%  data:2020年10月11日  author:飞蓬大将军 %%%%%%%%%%

%********************** 程序主体 ************%

function [demodata] = qpskdemod(idata,qdata,para,nd,ml)

demodata = zeros(para,ml*nd);
demodata((1:para),(1:ml:ml*nd-1)) = idata((1:para),(1:nd))>=0;
demodata((1:para),(2:ml:ml*nd)) = qdata((1:para),(1:nd))>=0;

end



