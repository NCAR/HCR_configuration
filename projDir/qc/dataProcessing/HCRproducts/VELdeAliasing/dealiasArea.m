function velUnfold=unfoldVel(velFolded,flag,elev)
% Unfold velocities

foldThresh=7.8311;

velUnfold=nan(size(velFolded));

velFolded(flag~=1)=nan;
velFolded(:,elev<0)=-velFolded(:,elev<0);

velGray=mat2gray(velFolded,[-foldThresh foldThresh]);

velHigh=velFolded>6.5;
end