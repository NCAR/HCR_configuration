function pid_comb= combine_pid_hcr_hsrl(pid_hcr,pid_hsrl,temp)



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

% Put HSRL SLW in everywhere

pid_comb(pid_hsrl==2)=2;
pid_comb(pid_hsrl==5)=5;

%pid_comb(pid_hsrl==6)=5; % HSRL SLW overrides HCR
%pid_comb(pid_hsrl==7)=3; % HSRL Drizzle overrides HCR


% Put HSRL ice in everywhere except HCR identified aggregates
%pid_hsrl(pid_comb==7)=0; % HCR Aggregates overrides HSRL ice
%pid_comb(pid_hsrl==5)=5;%pid_comb(pid_hsrl==5)=6; % HSRL ice overrides HCR
%pid_comb(pid_hsrl==4)=6; % HSRL ice overrides HCR
pid_comb(pid_hsrl==6)=6; % HSRL ice overrides HCR


% Put HSRL aerosol in everywhere except where HCR detects data
%pid_hsrl(isnan(pid_hcr)==0)=0;
%pid_comb((pid_hcr==1) | (pid_hcr==nan))=9;

%  pid_comb(pid_hsrl==4 & pid_comb==1)=9;
%  pid_comb(pid_hsrl==7 & pid_comb==1)=9;

pid_temp=ones(size(pid_hsrl));
pid_temp(pid_hsrl==4)=9;
pid_temb(pid_hsrl==7)=9;
pid_comb(pid_comb==1)=pid_temp(pid_comb==1);

% pid_comb(pid_hsrl==2 & pid_hcr > 8)=9;%
% pid_comb(pid_hsrl==3 & pid_hcr > 8)=9;%
end



