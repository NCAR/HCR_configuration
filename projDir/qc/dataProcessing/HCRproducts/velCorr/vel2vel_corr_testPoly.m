function velSmooth = vel2vel_corr_poly(VEL_ground,time,ptpIn,polyOrder)
% This function calculates VEL_CORR from VEL (it doesn't do the correction
% for the motion, i.e. VEL from VEL_RAW)
%
% Input variables:
% VEL: motion corrected velocity
% time
% linInd: linear index of ocean surface
%
% Output variables
% VEL_CORR
% velSmooth: polinomial fit curve

ptp=floor(ptpIn*10/2);

velSmooth=nan(length(VEL_ground),1);

gapNan=0;

for ii=(ptp+1):length(time)-ptp
    dataInd=ii-ptp:1:ii+ptp;
    grabInd=ptp+1;
    
    timeNum=datenum(time(dataInd));
    
    [polFit S Mu]= polyfit(timeNum,VEL_ground(dataInd),polyOrder);
    
    if any(isnan(polFit))
        gapNan=1;
    end
    
    if ~any(isnan(polFit)) & gapNan==1
        velSmooth(ii-ptp:ii+ptp)=polyval(polFit,timeNum,[],Mu);
        gapNan=0;
    elseif gapNan==0
        velSmooth(ii:ii+ptp)=polyval(polFit,timeNum(grabInd:grabInd+ptp),[],Mu);
    end
end
end

