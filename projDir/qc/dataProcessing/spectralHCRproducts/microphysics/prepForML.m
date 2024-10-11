function dataForML=prepForML(dataIn,vars)

% Set limits for scaling
lims.DBZ=[-50,20];
lims.VEL=[-18,18];
lims.WIDTH=[0,3];
lims.SKEW=[-3,3];
lims.KURT=[-8,8];
lims.EE_WIDTH=[0,15];
lims.L_SLOPE=[0,25];
lims.R_SLOPE=[-25,0];
lims.MELTING_LAYER=[9,21];

stackVars=[];

% Scale
for ii=1:length(vars)
    varThis=vars{ii};

    varThis=dataIn.(varThis);
    varThis(varThis<lims.(vars{ii})(1))=lims.(vars{ii})(1);
    varThis(varThis>lims.(vars{ii})(2))=lims.(vars{ii})(2);
    varScaled=1/(lims.(vars{ii})(2)-lims.(vars{ii})(1)).*(varThis-lims.(vars{ii})(1));

    % colmap=turbo(10);
    % surf(varScaled,'EdgeColor','none');
    % colormap(colmap)
    % view(2)
    % colorbar

    stackVars=cat(3,stackVars,varScaled);
end

dataForML=double(reshape(stackVars,size(dataIn.(vars{ii}),1)*size(dataIn.(vars{ii}),2),[],1));

end