function diffBoth=findFoldBoundaries(velFolded,testCrit)
diffHor=diff(velFolded,1,2);
diffHor=cat(2,nan(size(velFolded,1),1),diffHor);
diffVer=diff(velFolded,1,1);
diffVer=cat(1,nan(1,size(velFolded,2)),diffVer);

diffBoth=zeros(size(diffVer));

diffBoth(abs(diffHor)>testCrit)=1;
diffBoth(abs(diffVer)>testCrit)=1;

% Connect lines
diffBoth=bwmorph(diffBoth,'bridge');

% Remove small lines
diffBoth=bwareaopen(diffBoth,5);
end