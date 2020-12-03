function [pid_comb which_pid]= combine_pid_hcr_hsrl_clean(pid_hcr,pid_hsrl)

% pid_hsrl
%  nan   no signal
%  1   cloud liquid
%  2   Drizzle
%  3   Aerosol1 
%  4   SLW
%  5   Ice crystals
%  6   Aerosol2

% pid_hcr
%  nan no signal
%  1 cloud liquid
%  2 Drizzle
%  3 Rain
%  4 SLW
%  5 Ice crystals
%  6 Snow
%  7 wet snow/rimed ice

% pid_comb
%  nan no signal
%  1 cloud liquid
%  2 Drizzle
%  3 Rain
%  4 SLW
%  5 Ice crystals
%  6  Snow
%  7 wet snow/rimed ice
%  8 Aerosols


pid_comb=pid_hcr;
which_pid=zeros(size(pid_comb));

% Put HSRL SLW in everywhere
pid_comb(pid_hsrl==1)=1;
pid_comb(pid_hsrl==4)=4;
which_pid(pid_hsrl==1)=1;
which_pid(pid_hsrl==4)=1;

% This is something Vivek put in in the latest version but it is strange so
% I am not using it for now
%pid_comb(pid_hsrl==2 & isnan(pid_hcr)=8;
%pid_comb(pid_hsrl==1 & isnan(pid_hcr))=8;

% Put HSRL ice in everywhere except HCR identified aggregates
pid_comb(pid_hsrl==5)=5; % HSRL ice overrides HCR
which_pid(pid_hsrl==5)=1;

% Put HSRL aerosol in everywhere except where HCR detects data
pid_comb(pid_hsrl==3 & isnan(pid_comb))=8;
pid_comb(pid_hsrl==6 & isnan(pid_comb))=8;
which_pid(pid_comb==8)=1;
which_pid(pid_comb==8)=1;

end



