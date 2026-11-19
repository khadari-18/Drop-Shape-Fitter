%%% subfunction for derivative iterations in runge-kutta schemes %%%
function [dds_phi,dds_x,dds_y] = dds(phi,x,y,radius,bond)
    
%% the localized dervatives of spatial parameters %%

dds_phi = 2/radius - bond*y/radius^2 - sin(phi)/x;
% dds_phi = 2 - bond*y - sin(phi)/x;
dds_x = cos(phi);
dds_y = sin(phi);

%% commensts on issues @ apex with the following ways %%
% 
% ddx_phi = 2/cos(phi) - Bo*z/cos(phi) - tan(phi)/x; 
% ddx_z = tan(phi);    
% 
% ddz_phi = 2/sin(phi) - Bo*z/sin(phi) - 1/x; 
% ddz_x = cot(phi);

end