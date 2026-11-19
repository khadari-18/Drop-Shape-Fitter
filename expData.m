%%% outline extraction from image %%%
function [recon,coord_m,count,contour_m,Ro_m,px_by_m] = expData(im,needle_dia_m,aspect_ratio,tr)

% img = imread('PBSv204[9].JPG');
% needle_dia_m = 0.0018;
% tr = 0.06;

im = im2gray(im);                                                           % grayscale conv.
level = graythresh(im);                                                     % between 0 and 1

im = ~imbinarize(im,level);                                                 % binary conv.
im = imfill(im,"holes");

im = edge(im,"Canny",level);

col_half = round(size(im,2)/2);

im1 = im(:,1:col_half);
im2 = im(:,col_half+1:end);

[r1, c1] = find(im1==1);
[r2, c2] = find(im2==1);

[y1,i1] = sort(r1, "ascend");
[y2,i2] = sort(r2, "descend");

x1 = c1(i1);
x2 = c2(i2)+col_half;

recon = [x1,y1; x2,y2];

recon = smoothdata(recon, "gaussian", 32);                                  % smoothing contour shape from pixel induced zigzags 

%     imshow(im);
%     hold on
%     plot(recon(:,1), recon(:,2),'green');
    

%% pixel to world coordinate ratio %%
needle_dia_px = norm(recon(end,1)-recon(1,1));                              % unit conversion ratio      
px_by_m = needle_dia_px/needle_dia_m;                                       % the ratio
% px_by_m = (10^6)/9.772;

% trimming near needle data points at the top
count = size(recon,1);
recon = recon(round(tr*count/2):count-round(tr*count/2),:);                 % trim top-edge points, too close to needle
count = size(recon,1);

%% unit conversion and shifting origin to apex
coord_m = recon/px_by_m;                                                    % converted to metric

[apx_y, apx_ind] = max(coord_m(:,2));                                         
apx_x = coord_m(apx_ind,1);

coord_m(:,1) = coord_m(:,1)-apx_x;                                          % bringing apex on x = 0 line
coord_m(:,2) = abs(coord_m(:,2)-apx_y)/aspect_ratio;                        % capsize to upright pendant in the +ve y-axis (plus fix for the aspect ratio)


%% contour length and Ro guess in m 
contour_m = 0;                                                              % contour length storage init.

for i = 2:count
	contour_m = contour_m+norm(coord_m(i,:)-coord_m(i-1,:));
end

R_frac = 0.85;                                                              % 0.85 --> to initially consider bottom 15% of data for R-fit     
data_Ro = coord_m(round(R_frac*count/2)...
    :count-round(R_frac*count/2),:);
Ro_m = R_o(data_Ro(:,1), data_Ro(:,2));                                     % R of largest circle fitted @ apex 

%% drop volume calc. %%
% vol = 0;                                                                    % drop volume storage init.
% i = 0;                                                                      % reset counter
% while i ~= round(count/2)-1
%     vol = vol+0.5* ...
%         ((coord(end-i-1,1)-coord(i+2,1))^2+ ...
%         (coord(end-i,1)-coord(i+1,1))^2)* ...
%         (coord(i+2,2)-coord(i+1,2));
%     i = i+1;
% end
% 
% vol = 3.1416*vol;                                                           % in m^3
%
% imshow(im_org);
% hold on
% plot(recon(:,1), recon(:,2),'green');

end