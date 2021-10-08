function [rawDir interp] = modelDir(project,model,qc,qcVersion,freq)
% Find model directory
if strcmp(project,'socrates')
    baseDir='/scr/snow2/rsfdata/projects/socrates/';
elseif strcmp(project,'cset')
    baseDir='/scr/snow2/rsfdata/projects/cset/';
elseif strcmp(project,'otrec')
    baseDir='/scr/sleet2/rsfdata/projects/otrec/';
elseif strcmp(project,'spicule')
    baseDir='/scr/sleet2/rsfdata/projects/spicule/';
end

rawDir=[baseDir,'model/',model,'/'];
interp=[baseDir,'hcr/',qc,'/',model,'interp/',qcVersion,'/',freq,'/'];
end