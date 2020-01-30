%gui_plot_hcr  : 17:34 July 9th 2015


DBMVC=plot_ground_HCR(DBMVC,DBZ,altitude,range,time);
sur_sigma0=HCR_sur_cons(DBMVC,altitude,elevation);

grid

% Function var_ground
function var_ground=plot_ground_HCR(var,DBZ,altitude,range,time)

[m,n]=size(var);

% find range index of the ground echo using max dBZ
DBZ_tmp=DBZ;
DBZ_tmp(:,1:15)=0;
[~,ground_index]=max(DBZ_tmp,[],2);clear DBZ_tmp

% initialize var_ground
var_ground(1:m)=nan;

% populate var_ground
for i=1:m
    if abs(altitude(i)-range(ground_index(i))) < 0.5*altitude(i)
        var_ground(i)=var(i,ground_index(i));
    else
        var_ground(i)=nan;
    end
end

figure;plot(time,var_ground)
xlabel('Time, sec','FontSize',14);

end
% Function surface sigma0

    function sur_sigma0=HCR_sur_cons(dBm,alt,ele)
        % Peak transmit power=1 kW = 60 dBm
        Pr=10.^((dBm)*0.1); % milli watts
        pulse_width=0.256e-6; %  second; prf in Hz; avg_power in mWatts
        wave_len=0.00317797e-3; %transmit wavelength in km;
        
        RC=10^(-7.25);  % HCR radar contant
        c=3e8; % m/sec
        
        rc1=c*(pi^5)*pulse_width/(2.*(wave_len^4)*1e21);
        
        gaseous_att= 0; % two-way loss in db
        
        
        sur_sigma0=10*log10(RC*rc1*Pr.*(alt.^2)./sin(pi*ele./180))-gaseous_att;
        
        % The above is eq.(3) from Lihau Li et al. JTECH 2005.
        % In the denominator sine is used instead of cosine to be consistent with angles written in HCR data
        % See the adjustment to the angle when te plotting is done in the command
        % below
        figure;plot(ele+90,sur_sigma0)
        
        title('Sigma0 vs elevation angle','FontSize',18,'fontweight','b');
        xlabel('Elevation, degree','FontSize',16,'fontweight','b')
        ylabel('Normalized Sigma0, dB','FontSize',16,'fontweight','b')
    end