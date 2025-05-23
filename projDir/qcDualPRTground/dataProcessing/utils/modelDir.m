function [rawDir interp] = modelDir(project,model,qc,qcVersion,freq)
% Find model directory
if strcmp(project,'socrates')
    baseDir='/scr/snow2/rsfdata/projects/socrates/';
elseif strcmp(project,'cset')
    baseDir='/scr/snow2/rsfdata/projects/cset/';
elseif strcmp(project,'otrec')
    baseDir='/scr/sleet2/rsfdata/projects/otrec/';
elseif strcmp(project,'spicule')
    baseDir='/scr/sleet3/rsfdata/projects/spicule/';
elseif strcmp(project,'noreaster')
    baseDir='/scr/snow2/rsfdata/projects/noreaster/';
elseif strcmp(project,'meow')
    baseDir='/scr/virga1/rsfdata/projects/meow/';
end

rawDir=[baseDir,'model/',model,'/'];
if strcmp(freq,'10hz') | strcmp(freq,'2hz') | strcmp(freq,'10hz_combined')
    interp=[baseDir,'hcr/',qc,'/',model,'interp/',qcVersion,'/',freq,'/'];
elseif strcmp(freq,'combined')
    interp=[baseDir,'hcr/',qc,'/',model,'interp/',qcVersion,'/hcr_hsrl_merge/2hz/'];
end
end