function topoDirectory = topoDir(project)
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

topoDirectory=[baseDir,'topo/'];
end