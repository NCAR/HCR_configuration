function [pid_comb which_pid]= combine_pid_hcr_hsrl_clean(pid_hcr,pid_hsrl)

% pid_hsrl
%  1   no signal
%  2   cloud liquid
%  3   Drizzle
%  4   Aerosol1 
%  5   SLW
%  6   Ice crystals
%  7   Aerosol2

% pid_hcr
%  1 no signal
%  2 cloud liquid
%  3 Drizzle
%  4 Rain
%  5 SLW
%  6 Ice crystals
%  7 Snow
%  8 wet snow/rimed ice

% pid_comb
%  1 no signal
%  2 cloud liquid
%  3 Drizzle
%  4 Rain
%  5 SLW
%  6 Ice crystals
%  7  Snow
%  8 wet snow/rimed ice
%  9 Aerosols


pid_comb=pid_hcr;
which_pid=zeros(size(pid_comb));

% Put HSRL SLW in everywhere
pid_comb(pid_hsrl==2)=2;
pid_comb(pid_hsrl==5)=5;
which_pid(pid_hsrl==2)=1;
which_pid(pid_hsrl==5)=1;

% Put HSRL ice in everywhere except HCR identified aggregates
pid_comb(pid_hsrl==6)=6; % HSRL ice overrides HCR
which_pid(pid_hsrl==6)=1;

% Put HSRL aerosol in everywhere except where HCR detects data
pid_comb(pid_hsrl==4 & pid_comb==1)=9;
pid_comb(pid_hsrl==7 & pid_comb==1)=9;
which_pid(pid_comb==9)=1;
which_pid(pid_comb==9)=1;

% pid_temp=ones(size(pid_hsrl));
% pid_temp(pid_hsrl==4)=9;
% pid_temb(pid_hsrl==7)=9;
% pid_comb(pid_comb==1)=pid_temp(pid_comb==1);

end



