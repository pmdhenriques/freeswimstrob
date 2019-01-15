function Orientation_C = fish_orientation(Orientation,Centroids,BBoxs)
% Computes the real orientation of a skewd body in degrees 
% (like that of a fish) from the Orientation, Bounding box
% and centroid positions of the object
Orientation_C = NaN(length(Orientation),1);

boxW = BBoxs(:,3); % BB Width
boxH = BBoxs(:,4); % BB Height
boxX = BBoxs(:,1); % BB X
boxY = BBoxs(:,2); % BB Y

X = bsxfun(@minus,bsxfun(@plus,boxX,(boxW./2)),Centroids(:,1));  % Distance of centroid xpos from the center of the BB
Y = bsxfun(@minus,bsxfun(@plus,boxY,(boxH./2)),Centroids(:,2));  % Distance of centroid ypos from the center of the BB

% Define in which quadrant of the BB the centroid is and add the propper
% angle

ix1 = X < 0 & Y >= 0;
ix2 = X < 0 & Y < 0;
ix3 = X >= 0 & Y < 0;
ix4 = X >= 0 & Y >= 0;

Orientation_C(ix1) = bsxfun(@minus,90,Orientation(ix1));
Orientation_C(ix2) = bsxfun(@plus,abs(Orientation(ix2)),90);
Orientation_C(ix3) = bsxfun(@plus,bsxfun(@minus,Orientation(ix3),90),180);
Orientation_C(ix4) = bsxfun(@plus,abs(Orientation(ix4)),270);

end