function filename=nctest2(theNetCDFFile)

% Open file.
ncid = netcdf.open(theNetCDFFile,'NC_NOWRITE');

% Get number of dimensions, number of variables, number of global 
% attributes, and the identity of the unlimited dimension
[~,nvars,~,~] = netcdf.inq(ncid);

% Loop over number of variables
for i=1:nvars
    % Get the name, datatype, dimensions IDs, and the number of attributes 
    % of the variable.
    [varname, xtype, ~, ~] = netcdf.inqVar(ncid,i-1);
    data=ncread(theNetCDFFile,varname);
    assignin('caller',varname,data');
end

    % Get variable ID of the variable, given its name.
%      varid = netcdf.inqVarID(ncid,varname);
%     
%     % Get the value of the variable, given its ID.
%     % First convert single precision (xtype=5) to double otherwise 
%     % read it as it exists.
%     if xtype==5 
%         data = netcdf.getVar(ncid,varid,'double');
%     else 
%         data = netcdf.getVar(ncid,varid);  
%     end
%     
%     % apply scale factor and offset if they are defined and convert to
%     % double.
%     scale_factor=0.0;
%     add_offset=0.0;
%     missing_value=-9e10;
%     for j=1:varAtts
%         attname=netcdf.inqAttName(ncid,varid,j-1);
%         if strcmpi(attname,'scale_factor')
%             scale_factor = netcdf.getAtt(ncid,varid,attname,'double');
%         end
%         if strcmpi(attname,'add_offset')
%             add_offset = netcdf.getAtt(ncid,varid,attname,'double');
%         end
%         if strcmpi(attname,'_FillValue')
%             missing_value = netcdf.getAtt(ncid,varid,attname,'double');
%         end
%     end
% 
%     if scale_factor ~= 0.0
%         data=double(data);
%         % don't scale missing values
%         data(data<missing_value+0.10)=nan;
%         data=data*scale_factor + add_offset;
%     end
%     % Assign data to the variable varname
%     assignin('caller',varname,data');
% end