%%%%%%%%%%%%%%%%%%%%% OFDM仿真 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% fade.m  %%%%%%%%%
%%%%%%%%%  data:2020年10月17日  author:飞蓬大将军 %%%%%%%%%%

% idata : input Ich data
% qdata : input Qch data
% iout : output Ich data
% qout : output Qch data
% ramp : Amplitude contaminated by fading
% rcos : Cosine value contaminated by fading
% rsin : Cosine value contaminated by fading
% nsamp : Number of samples to be simulated
% tstp : Minimum time resolution
% fd : maximum doppler frequency
% no : number of waves in order to generate fading
% counter : fading counter
% flat : flat fading or not

function [iout,qout,ramp,rcos,rsin] = fade(idata,qdata,nsamp,tstp,fd,no,counter,flat)
if fd~=0.0
    ac0 = sqrt(1.0./(2.0.*(no+1)));
    % power normalized constant(ich)
    as0 = sqrt(1.0./(2.0.*no));
    % power normalized constant(qch)
    ic0 = counter;
    % fading counter
    pai = 3.14159265;
    wm = 2.0.*pai.*fd;
    n = 4.*no + 2;
    ts = tstp;
    wmts = wm.*ts;
    paino = pai./no;
    
    xc = zeros(1,nsamp);
    xs = zeros(1,nsamp);
    ic = [1:nsamp] +ic0;
    
    for nn =1:no
       cwn = cos(cos(2.0.*pai.*nn./n).*ic.*wmts); %cwn是什么？
       xc = xc + cos(paino.*nn).*cwn;
       xs = xs + cos(paino.*nn).*cwn;
    end
    cwmt = sqrt(2.0).*cos(ic.*wmts);
    xc = (2.0.*xc + cwmt).*ac0; %xc是什么？
    xs = 2.0.*xs.*as0;
    
    ramp = sqrt(xc.^2+xs.^2);
    rcos  = xc./ramp;
    rsin = xs./ramp;
    
    if flat == 1
        iout = sqrt(xc.^2+xs.^2).*idata(1:nsamp);
        qout = sqrt(xc.^2+xs.^2).*qdata(1:nsamp);
    else
        iout = xc.*idata(1:nsamp) - xs.*qdata(1:nsamp); 
        qout = xs.*idata(1:nsamp) + xc.*qdata(1:nsamp);
    end
else
    iout = idata;
    qout = qdata;
        
        
end
    



end