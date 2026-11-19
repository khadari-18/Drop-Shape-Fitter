%%% rk5 for Y-L equation to generate theoretical drop-outline %%%
function [theo_outline] = theoryData(contour_m,radius,bond) 

num = 4999;
s_step = contour_m/(num-1);

%% data storage and initial conditions 
coord_m = zeros(num,3);                                                     % phi (in radians everywhere), x, and y storage 

%% Runge-Kutta loops 
for i=(num+1)/2:num-1

	if i==(num+1)/2
        J1 = 1/radius;                                                      % dphi/ds initial value
        K1 = 1;                                                             % dx/ds initial value
        L1 = 0;                                                             % dy/ds initial value
	elseif i>1
        [J1, K1, L1] = ...
            dds(coord_m(i,1),coord_m(i,2),coord_m(i,3),radius,bond);
	end

%% RK4 steps
    [J2, K2, L2] = dds(coord_m(i,1)+J1*s_step/2,coord_m(i,2)+K1*s_step/2, ...
        coord_m(i,3)+L1*s_step/2,radius,bond);
    [J3, K3, L3] = dds(coord_m(i,1)+J2*s_step/2,coord_m(i,2)+K2*s_step/2, ...
        coord_m(i,3)+L2*s_step/2,radius,bond);
    [J4, K4, L4] = dds(coord_m(i,1)+J3*s_step,coord_m(i,2)+K3*s_step, ...
        coord_m(i,3)+L3*s_step,radius,bond);

    coord_m(i+1,1) = coord_m(i,1) + (J1+2*J2+2*J3+J4)*s_step/6; 
    coord_m(i+1,2) = coord_m(i,2) + (K1+2*K2+2*K3+K4)*s_step/6; 
    coord_m(i+1,3) = coord_m(i,3) + (L1+2*L2+2*L3+L4)*s_step/6;    

% %% RK5 steps
%     [J2, K2, L2] = dds(coord(i,1)+J1*s_step/4,coord(i,2) ...
%         +K1*s_step/4, coord(i,3)+L1*s_step/4, Bo, Ro);
%     [J3, K3, L3] = dds(coord(i,1)+J1*s_step/8+J2*s_step/8,coord(i,2) ...
%         +K1*s_step/8+K2*s_step/8, coord(i,3)+L1*s_step/8 ...
%         +L2*s_step/8, Bo, Ro);
%     [J4, K4, L4] = dds(coord(i,1)-J2*s_step/2+J3*s_step, ...
%         coord(i,2)-K2*s_step/2+K3*s_step, ...
%         coord(i,3)-L2*s_step/2+L3*s_step, Bo, Ro);
%     [J5, K5, L5] = dds(coord(i,1)+J1*s_step*3/16+J4*s_step*9/16, ...
%         coord(i,2)+K1*s_step*3/16+K4*s_step*9/16, ...
%         coord(i,3)+L1*s_step*3/16+L4*s_step*9/16, Bo, Ro);
%     [J6, K6, L6] = dds(coord(i,1)-J1*s_step*3/7+J2*s_step*2/7 ...
%         +J3*s_step*12/7-J4*s_step*12/7+J5*s_step*8/7, ...
%         coord(i,2)-K1*s_step*3/7+K2*s_step*2/7+K3*s_step*12/7 ...
%         -K4*s_step*12/7+K5*s_step*8/7, coord(i,3)-L1*s_step*3/7 ...
%         +L2*s_step*2/7+L3*s_step*12/7-L4*s_step*12/7 ...
%         +L5*s_step*8/7, Bo, Ro);
% 
%     coord(i+1,1) = coord(i,1) + (7*J1+32*J3+12*J4+32*J5+7*J6)*s_step/90; 
%     coord(i+1,2) = coord(i,2) + (7*K1+32*K3+12*K4+32*K5+7*K6)*s_step/90; 
%     coord(i+1,3) = coord(i,3) + (7*L1+32*L3+12*L4+32*L5+7*L6)*s_step/90;

end


%% creating full shape and processing data for co-plotting with image 
theo_outline = zeros(num,2);

% creating the full shape first
theo_outline(:,1) = ...
    [-flip(coord_m((num+1)/2:num,2)); coord_m(1+(num+1)/2:num,2)];          % x-coordinates                  
theo_outline(:,2) = ...
    [flip(coord_m((num+1)/2:num,3)); coord_m(1+(num+1)/2:num,3)];           % y-coordinates                  


%% viewing finally processed data 
% subplot(1,2,1)
% plot(theo_outline(:,1),theo_outline(:,2));

% subplot(1,2,2)
% plot(coord(2,:),coord(3,:));   

end
