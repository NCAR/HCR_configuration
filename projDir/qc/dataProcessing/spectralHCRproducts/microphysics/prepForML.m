function dataForML=prepForML(dataIn,vars)

% Set limits for scaling
lims.DBZ_MASKED=[-50,20];
lims.VEL_MASKED=[-18,18];
lims.WIDTH_SPEC=[0,3];
lims.SKEWNESS=[-3,3];
lims.KURTOSIS=[-8,8];
lims.EDGE_EDGE_WIDTH=[0,15];
lims.LEFT_SLOPE=[0,25];
lims.RIGHT_SLOPE=[-25,0];
lims.MELTING_LAYER=[9,21];
lims.TEMP=[-80,50];

stackVars=[];

% Scale
for ii=1:length(vars)
    varThis=vars{ii};

    varThis=dataIn.(varThis);
    varThis(varThis<lims.(vars{ii})(1))=lims.(vars{ii})(1);
    varThis(varThis>lims.(vars{ii})(2))=lims.(vars{ii})(2);
    varScaled=1/(lims.(vars{ii})(2)-lims.(vars{ii})(1)).*(varThis-lims.(vars{ii})(1));

    stackVars=cat(3,stackVars,varScaled);
end

dataForML=double(reshape(stackVars,size(dataIn.(vars{ii}),1)*size(dataIn.(vars{ii}),2),[],1));

end