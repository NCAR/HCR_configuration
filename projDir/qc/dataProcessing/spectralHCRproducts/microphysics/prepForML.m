function [dataForML,lims]=prepForML(dataIn,vars)

% Set limits for scaling
lims.DBZ_MASKED=[-50,20];
lims.VEL_MASKED=[-10,10];
lims.WIDTH_SPEC=[0,2.5];
lims.SKEWNESS=[-2,2];
lims.KURTOSIS=[-6,6];
lims.EDGE_EDGE_WIDTH=[0,15];
lims.LEFT_SLOPE=[0,20];
lims.RIGHT_SLOPE=[-20,0];
lims.LEFT_EDGE_VEL=[-12,12];
lims.RIGHT_EDGE_VEL=[-12,12];
lims.MELTING_LAYER=[9,21];
lims.TEMP=[-60,40];

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