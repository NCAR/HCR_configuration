function indir = HCRdirWFH(project,quality,freq)
% Find HCR data directory
%% SOCRATES
if strcmp(project,'socrates')
    if strcmp(quality,'qc2')
        if strcmp(freq,'2hz') | strcmp(freq,'10hz') | strcmp(freq,'100hz')
            indir=['/run/media/romatsch/RSF0006/rsf/hcr/socrates/',freq,'/'];
        elseif strcmp(freq,'combined')
            indir='/run/media/romatsch/RSF0006/rsf/combined_hcr_hsrl/socrates/';
        else
            disp('Frequency not available.')
            return
        end
    else
        disp('Only qc2 data available.')
        return
    end
    
    %% CSET
elseif strcmp(project,'cset')
    if strcmp(quality,'qc2')
        if strcmp(freq,'2hz') | strcmp(freq,'10hz')
            indir=['/run/media/romatsch/RSF0006/rsf/hcr/cset/',freq,'/'];
        elseif strcmp(freq,'combined')
            indir='/run/media/romatsch/RSF0006/rsf/combined_hcr_hsrl/cset/';
        else
            disp('Frequency not available.')
            return
        end
    else
        disp('Only qc2 data available.')
        return
    end
    
    %% OTREC
elseif strcmp(project,'otrec')
    if strcmp(quality,'qc2')
        if strcmp(freq,'10hz')
            indir=['/run/media/romatsch/RSF0006/rsf/hcr/otrec/',freq,'/'];
        else
            disp('Frequency not available.')
            return
        end
    else
        disp('Only qc2 data available.')
        return
    end
end
end

