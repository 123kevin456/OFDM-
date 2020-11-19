%%%%%%%%%%%%%%%%%%%%% �ӳٺ��� %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% delay.m  %%%%%%%%%
%%%%%%%%%  data:2020��10��19��  author:����󽫾� %%%%%%%%%%


%%%%%����˵��
%%%�ӳٺ���

%%%%    ���滷��
%����汾��MATLAB R2019a
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