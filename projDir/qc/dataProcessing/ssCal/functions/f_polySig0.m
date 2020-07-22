function PLT = f_polySig0(PLT)
% Make poly fit to sig0 model and data

elevPoly=0:0.5:20;

elevData=cat(2,PLT.elev,PLT.sig0measured);
elevData(any(isnan(elevData), 2), :) = [];
polyData=polyfit(elevData(:,1),elevData(:,2),3);
PLT.polySig0Meas = polyval(polyData,elevPoly);

if min(isnan(PLT.sig0model(:,2)))==0
    elevModel=cat(2,PLT.elev,PLT.sig0model(:,2));
    elevModel(any(isnan(elevModel), 2), :) = [];
    polyModel=polyfit(elevModel(:,1),elevModel(:,2),3);
    polySig0Model1 = polyval(polyModel,elevPoly);
else
    polySig0Model1=nan(size(elevPoly));
end

if min(isnan(PLT.sig0model(:,5)))==0
    elevModel=cat(2,PLT.elev,PLT.sig0model(:,5));
    elevModel(any(isnan(elevModel), 2), :) = [];
    polyModel=polyfit(elevModel(:,1),elevModel(:,2),3);
    polySig0Model2 = polyval(polyModel,elevPoly);
else
    polySig0Model2=nan(size(elevPoly));
end

if min(isnan(PLT.sig0model(:,8)))==0
    elevModel=cat(2,PLT.elev,PLT.sig0model(:,8));
    elevModel(any(isnan(elevModel), 2), :) = [];
    polyModel=polyfit(elevModel(:,1),elevModel(:,2),3);
    polySig0Model3 = polyval(polyModel,elevPoly);
else
    polySig0Model3=nan(size(elevPoly));
end

PLT.polySig0Model=cat(1,polySig0Model1,polySig0Model2,polySig0Model3);

PLT.polyElev=elevPoly;

end

