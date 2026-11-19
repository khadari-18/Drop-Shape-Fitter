%% Calibration for the vert/hor pixel aspect ratio %%
sphere = imread('Calib 5.png');
sphere = im2gray(sphere);
level = graythresh(sphere);

sphere = ~imbinarize(sphere,level);

bound_vert = round(0.42*size(sphere,1)):round(0.58*size(sphere,1));
bound_hor = round(0.42*size(sphere,2)):round(0.58*size(sphere,2));

% bound_vert = round(0.46*size(sphere,1)):round(0.48*size(sphere,1));
% bound_hor = round(0.49*size(sphere,2)):round(0.51*size(sphere,2));

strip_vert = sphere(:,bound_vert);                                          
strip_hor = sphere(bound_hor,:);                                           

edge_vert = edge(strip_vert,"Canny",level);
edge_hor = edge(strip_hor,"Canny",level);

[coord_vert,~] = find(edge_vert==1);
[~,coord_hor] = find(edge_hor==1);

dia_px_vert = max(coord_vert)-min(coord_vert);
dia_px_hor = max(coord_hor)-min(coord_hor);

Aspect_ratio = dia_px_vert/dia_px_hor;


