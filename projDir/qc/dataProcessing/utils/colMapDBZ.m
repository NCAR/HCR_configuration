function [fig] = colMapDBZ(fig)
%Change dbz color map

fld=fig.CData;

col_def1 = nan(size(fld));
col_def2 = nan(size(fld));
col_def3 = nan(size(fld));

limits=[-inf (-43:3:23) inf];
color_map=colormap(dbz_default);

for ii=1:size(color_map,1)
    col_ind=find(fld>limits(ii) & fld<=limits(ii+1));
    col_def1(col_ind)=color_map(ii,1);
    col_def2(col_ind)=color_map(ii,2);
    col_def3(col_ind)=color_map(ii,3);
end
if ~isequal(size(col_def1),(size(fld)))
    col_def=cat(3,col_def1',col_def2',col_def3');
else
    col_def=cat(3,col_def1,col_def2,col_def3);
end
fig.CData=col_def;

hcb=colorbar;
set(get(hcb,'Title'),'String','dBZ');
colormap(gca,color_map);
caxis([0 size(color_map,1)]);
caxis_yticks=(1:1:size(color_map,1)-1);
caxis_ytick_labels=num2str(limits(2:end-1)');
while length(caxis_yticks)>16
    caxis_yticks=caxis_yticks(1:2:end);
    caxis_ytick_labels=caxis_ytick_labels((1:2:end),:);
end
set(hcb,'ytick',caxis_yticks);
set(hcb,'YTickLabel',caxis_ytick_labels);
end

