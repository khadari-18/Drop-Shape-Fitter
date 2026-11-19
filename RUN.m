%% Input directory, needle diameter and image-top-trim inputs %%

host_path = 'D:\WORK\Processed data\Large Airway Mucus Characterization\Surface\DSA processed\Tween 20 DI\New Zoom Lens\01102025 (20 Gauge)\Conc. Sweep\0.001mM tseries\';
ndl_dia_m = 0.0008905;                                                      % input needle dia in meters here
aspect_ratio = 1;                                                           % vert_px / hor_px
trim_top = 0.15;                                                            % optional: between 1% to 10%

mkdir (host_path + "Results");                                              % creating Results folder

%% Loop to read the images named in the textfile in series

Impath = dir(append(host_path, '*.jpg'));                                   % contains only the file names, not paths included 
output_mat = cell(length(Impath)+1,4);                                      % to store name, st, area, vol 
output_mat(1,:) = {'name', 'sigma', 's_area', 'vol'};
k=1;
while k <= length(Impath)
    fullFileName = fullfile(host_path, Impath(k).name);

    % Imname = fgetl(ifid);                                                 % reading image names from text file
    % 
    % if Imname == -1 
    %     break;
    % end

    im = imread(fullFileName);                                              % input image title here

    %% experimental data extraction block %%
    [exp_px,exp_m,num,contour_m,Ro_m,px_by_m] = expData(im,ndl_dia_m,aspect_ratio,trim_top); % in metric, data pts always odd size returned 

    %% fitting block %%

    bethe = [0,0,0,Ro_m,0.3];                                               % optimization parameters = [X0,Y0,omega,Ro,Bo]

    % rotate experimental data for misalignment with gravity % 
    Ex_rot_handl = @(X0,Y0,omega) [(exp_m(:,1)-X0)*cos(omega)+(exp_m(:,2)-Y0)*sin(omega),(exp_m(:,1)-X0)*sin(omega)+(exp_m(:,2)-Y0)*cos(omega)];                                  

    % function handle defining objective that is minimized %
    resn = @(bethe) (sum(sqdistSearch(theoryData(contour_m,bethe(4),bethe(5)),feval(Ex_rot_handl, bethe(1),bethe(2),bethe(3))))); % sum of squared distances/residues

    % minimization function calling with specified options %
    options = optimset('TolFun',1e-10,'TolX',1e-10);
    [bethe,res] = fminsearch(resn,bethe,options);                              

    %% surface tension calculation w/ the fitted parameter values %%
    sigma = 9810000*(bethe(4)^2)/bethe(5);                                  % surface tension in mN/m unit


    %% plotting block %%

    % rotated experimental data %
    exp_m = feval(Ex_rot_handl,bethe(1),bethe(2),bethe(3)); 


    % drop surface area %
    s_area = 0;
    vol = 0;

    for i = round(num/2):num-1
	    s_area = s_area + exp_m(i,1)*norm(exp_m(i,:)-exp_m(i+1,:));         % area/2*pi = integration of x*ds (ds = sqrt(dx^2+dz^2))
        vol = vol + exp_m(i,1)^2 * norm(exp_m(i,:)-exp_m(i+1,:));   
    end

    s_area = s_area*2*3.1416;
    vol = vol*3.1416;

    % final theoretical outline %
    fitted_m = theoryData(contour_m,bethe(4),bethe(5));

    fig = figure('Visible','off');
    % overlayed experimental data plot %
    subplot(6,4,[1 2 5 6 9 10 13 14 17 18 21 22])
    imshow(im);
    hold on
    plot(exp_px(:,1), exp_px(:,2),'green');

    % fit check plot %
    subplot(6,4,[7 8 11 12 15 16 19 20])
    plot(exp_m(:,1),exp_m(:,2),'DisplayName','exp shape','Color','blue');
    hold on
    plot(fitted_m(:,1),fitted_m(:,2),'DisplayName','fitted shape','Color','red');
    legend('Location','north');

    saveas(fig, host_path + "Results\" + Impath(k).name, 'jpg');            % save the fitted image
    
    close(fig);

    output_mat(k+1,:) = {Impath(k).name, sigma, s_area, vol};
    k = k+1; 
end                                                                         % while end

    
xcelpath = fullfile(host_path, 'Results', 'results.xlsx');
writecell(output_mat,xcelpath);

