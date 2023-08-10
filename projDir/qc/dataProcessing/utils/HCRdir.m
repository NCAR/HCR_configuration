function indir = HCRdir(project,qc,qcVersion,freq)
% Find HCR data directory
indir=[];
%% SOCRATES
if strcmp(project,'socrates')
    if strcmp(qc,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/',freq,'/'];
        else
            disp('No 2hz data in field data.');
            return
        end
    elseif strcmp(qc,'qc1')
        if strcmp(freq,'100hz')
            indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/moments/100hz/';
        elseif strcmp(freq,'2hz')  | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/',freq,'/'];
        end
    elseif strcmp(qc,'qc2')
        if strcmp(freq,'2hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/',qc,'/cfradial/',qcVersion,'/',freq,'/'];
        elseif strcmp(freq,'2hzMerged')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/',qc,'/cfradial/hcr_hsrl_merge/',qcVersion,'/2hz/'];
        else
            disp('The requested data does not exist.');
            return
        end
    elseif strcmp(qc,'qc3')
        if strcmp(freq,'10hz') | strcmp(freq,'2hz')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/',qc,'/cfradial/',qcVersion,'_full/',freq,'/'];
        elseif strcmp(freq,'combined')
            indir=['/scr/snow2/rsfdata/projects/socrates/hcr/',qc,'/cfradial/hcr_hsrl_merge/',qcVersion,'_full/2hz/'];
        else
            disp('The requested data does not exist.');
            return
        end
    elseif strcmp(qc,'ts')
        indir='/scr/snow2/rsfdata/projects/socrates/hcr/time_series_netcdf/';
    end
    
    %% CSET
elseif strcmp(project,'cset')
    if strcmp(qc,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/',freq,'/'];
        else
            disp('No 2hz data in field data.');
            return
        end
    elseif strcmp(qc,'qc1')
        if strcmp(freq,'100hz')
            indir='/scr/snow2/rsfdata/projects/cset/hcr/qc/cfradial/moments/100hz/';
        elseif strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/qc/cfradial/velcorr/',freq,'/'];
        else
            disp('No 2hz data in cset data.');
            return
        end
    elseif strcmp(qc,'qc2')
        if strcmp(freq,'2hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/',qc,'/cfradial/',qcVersion,'/',freq,'/'];
        elseif strcmp(freq,'2hzMerged')
            indir='/scr/snow2/rsfdata/projects/cset/hcr/qc2/cfradial/hcr_hsrl_merge/',qcVersion,'_full/2hz/';
        else
            disp('The requested data does not exist.');
            return
        end
    elseif strcmp(qc,'qc3')
        if strcmp(freq,'10hz') | strcmp(freq,'2hz')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/',qc,'/cfradial/',qcVersion,'_full/',freq,'/'];
        elseif strcmp(freq,'combined')
            indir=['/scr/snow2/rsfdata/projects/cset/hcr/',qc,'/cfradial/hcr_hsrl_merge/',qcVersion,'_full/2hz/'];
        else
            disp('The requested data does not exist.');
            return
        end
    end
    
    %% ARISTO
elseif strcmp(project,'aristo')
    if strcmp(qc,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/',freq,'/'];
        else
            disp('No 2hz data in aristo data.');
            return
        end
    elseif strcmp(qc,'qc1') | strcmp(qc,'qc2')
        disp('There are no qc1 or qc2 data for aristo.')
        return
    end
    
    %% OTREC
elseif strcmp(project,'otrec')
    if strcmp(qc,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
           % indir=['/scr/snow1/rsfdata/projects/otrec/hcr/cfradial/moments/',freq,'/'];
           disp('Got lost.');
        end
    elseif strcmp(qc,'qc1')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            % indir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc1/cfradial/final/',freq,'/'];
            disp('Got lost.');
        else
            disp('No 2hz data in qc1 data.');
            return
        end
    elseif strcmp(qc,'qc2')
        if strcmp(freq,'10hz') & (strcmp(qcVersion,'v2.0') | strcmp(qcVersion,'v2.2'))
            indir=['/scr/sleet2/rsfdata/projects/otrec/hcr/',qc,'/cfradial/',qcVersion,'/',freq,'/'];
        else
            disp('The requested data does not exist.');
            return
        end
    elseif strcmp(qc,'qc3')
        if strcmp(freq,'10hz')
            indir=['/scr/sleet2/rsfdata/projects/otrec/hcr/',qc,'/cfradial/',qcVersion,'_full/',freq,'/'];
        else
            disp('The requested data does not exist.');
            return
        end
    elseif strcmp(qc,'ts')
        indir='/scr/sleet2/rsfdata/projects/otrec/hcr/time_series_netcdf/';
    end

    %% SPICULE
elseif strcmp(project,'spicule')
    if strcmp(qc,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/sleet3/rsfdata/projects/spicule/hcr/cfradial/moments/',freq,'/'];
        end
    elseif strcmp(qc,'qc0')
        if strcmp(freq,'10hz')
            indir=['/scr/sleet3/rsfdata/projects/spicule/hcr/',qc,'/cfradial/',qcVersion,'/',freq,'/'];
        end
    elseif strcmp(qc,'qc1')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/sleet3/rsfdata/projects/spicule/hcr/qc1/cfradial/',qcVersion,'_full/',freq,'/'];
        end
    elseif strcmp(qc,'ts')
        indir='/scr/sleet3/rsfdata/projects/spicule/hcr/time_series_netcdf/';
    end

    %% NOREASTER
elseif strcmp(project,'noreaster')
    if strcmp(qc,'field')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/noreaster/hcr/cfradial/moments/orig/',freq,'/'];
        end
%     elseif strcmp(qc,'qc0')
%         if strcmp(freq,'10hz')
%             indir=['/scr/sleet2/rsfdata/projects/spicule/hcr/',qc,'/cfradial/',qcVersion,'/',freq,'/'];
%         end
    elseif strcmp(qc,'qc1')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/noreaster/hcr/qc/cfradial/',qcVersion,'_full/',freq,'/'];
        end
    elseif strcmp(qc,'qc2')
        if strcmp(freq,'100hz') | strcmp(freq,'10hz')
            indir=['/scr/snow2/rsfdata/projects/noreaster/hcr/qc2/cfradial/',qcVersion,'_full/',freq,'/'];
        end
    end
end

end

