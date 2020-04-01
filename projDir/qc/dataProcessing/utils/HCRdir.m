function indir = HCRdir(project,quality,freq)
% Find HCR data directory
%% SOCRATES
if strcmp(project,'socrates')
    if strcmp(quality,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/',freq,'/'];
        else
            disp('No 2hz data in field data.');
            return
        end
    elseif strcmp(quality,'qc1')
        if strcmp(freq,'100hz')
            indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/moments/100hz/';
        elseif strcmp(freq,'2hz')  | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/',freq,'/'];
        end
    elseif strcmp(quality,'qc2')
        if strcmp(freq,'2hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/final/',freq,'/'];
        else
            disp('No 2hz data in qc2 data.');
            return
        end
    end
    
    %% CSET
elseif strcmp(project,'cset')
    if strcmp(quality,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/',freq,'/'];
        else
            disp('No 2hz data in field data.');
            return
        end
    elseif strcmp(quality,'qc1')
        if strcmp(freq,'100hz')
            indir='/scr/snow2/rsfdata/projects/cset/hcr/qc/cfradial/moments/100hz/';
        elseif strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/qc/cfradial/velcorr/',freq,'/'];
        else
            disp('No 2hz data in cset data.');
            return
        end
    elseif strcmp(quality,'qc2')
        if strcmp(freq,'2hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/qc2/cfradial/final/',freq,'/'];
        end
    end
    
    %% ARISTO
elseif strcmp(project,'aristo')
    if strcmp(quality,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/',freq,'/'];
        else
            disp('No 2hz data in aristo data.');
            return
        end
    elseif strcmp(quality,'qc1') | strcmp(quality,'qc2')
        disp('There are no qc1 or qc2 data for aristo.')
        return
    end
    
    %% OTREC
elseif strcmp(project,'otrec')
    if strcmp(quality,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow1/rsfdata/projects/otrec/hcr/cfradial/moments/',freq,'/'];
        end
    elseif strcmp(quality,'qc1')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc1/cfradial/final/',freq,'/'];
        else
            disp('No 2hz data in qc1 data.');
            return
        end
    elseif strcmp(quality,'qc2')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final/',freq,'/'];
        else
            disp('No 2hz data in qc1 data.');
            return
        end
    end
end
end

