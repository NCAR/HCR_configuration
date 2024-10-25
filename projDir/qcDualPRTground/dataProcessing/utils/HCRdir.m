function indir = HCRdir(project,qc,qcVersion,freq)
% Find HCR data directory
indir=[];
%% MEOW
if strcmp(project,'meow')
    if strcmp(qc,'ts')
        indir='/scr/virga1/rsfdata/projects/meow/hcr/time_series/wband/save/';
    elseif strcmp(qc,'qc0')
        indir=['/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/moments/',freq,'/'];
    else
        indir=['/scr/virga1/rsfdata/projects/meow/hcr/',qc,'/cfradial/',qcVersion,'_full/',freq,'/'];
    end
end

end

