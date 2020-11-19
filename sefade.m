%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% sefade.m  %%%%%%%%%
%%%%%%%%%  data:2020年10月16日  author:飞蓬大将军 %%%%%%%%%%


%%%%%程序说明
%%%频率选择性衰落

%%%%    仿真环境
%软件版本：MATLAB R2019a
% itau : Delay time for each multipath fading
% dlvl : Attenuation level for each multipath fading
% th : Initialized phase for each multipath fading
% n0 : Number of waves in order to generate each
% multipath fading
% itn : Fading counter for each multipath fading
% n1 : Number of summation for direct and delayed
% waves
% nsamp : Total number of symbols
% tstp : Minimum time resolution
% fd : Maximum doppler frequency
% flat flat fading or not
% (1-flat (only amplitude is fluctuated),0-normal(phase
% and amplitude are fluctuated))



function [iout,qout,ramp,rcos,rsin] = sefade(idata,qdata,itau,dlvl,th,n0,itn,nl,nsamp,tstp,fd,flat)


iout = zeros(1,nsamp);
qout = zeros(1,nsamp);
total_attn = sum(10.^(-1.0.*dlvl./10.0));

for k =1:nl
    atts = 10.^(-0.05.*dlvl(k));
    if dlvl(k) == 40.0
        atts = 0.0;
    end
    theta = th(k).*pi./180.0;
    [itmp,qtmp] = delay(idata,qdata,nsamp,itau(k));
    [itmp3,qtmp3,ramp,rcos,rsin] = fade(itmp,qtmp,nsamp,tstp,fd,n0(k),itn(k),flat);
    iout = iout + atts.*itmp3./sqrt(total_attn);
    qout = qout + atts.*qtmp3./sqrt(total_attn);
    
end % end for k


end