%%%%%%%%%%%%%%%%%%%%% 延迟函数 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% delay.m  %%%%%%%%%
%%%%%%%%%  data:2020年10月19日  author:飞蓬大将军 %%%%%%%%%%


%%%%%程序说明
%%%延迟函数

%%%%    仿真环境
%软件版本：MATLAB R2019a
function [iout,qout] = delay(idata,qdata,nsamp,idel)
iout = zeros(1,nsamp);
qout = zeros(1,nsamp);
if idel~=0
    iout(1:idel) = zeros(1,idel);
    qout(1:idel) = zeros(1,idel);
end
iout(idel+1:nsamp) = idata(1:nsamp-idel);
qout(idel+1:nsamp) = qdata(1:nsamp-idel);
   

end