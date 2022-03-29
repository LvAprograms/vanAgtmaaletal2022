%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to calculate slab pull force and %
% drag below overriding plate, as an      %
% extension of Ylona's vis_hdf5.m script  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luuk van Agtmaal, July 2020

clear all;
clfall = 2;
if clfall == 0
    % keep figure
elseif clfall == 1
    clf;
elseif clfall == 2
    fh=findall(0,'type','figure');
    for i=1:length(fh)
        clf(fh(i));
        close(fh(i));
    end
end

% -- PLOTTING SPECIFICATIONS --

% --- Which files ---
expname     = '/media/luuk/My Book/MSc_paper/Data/FI_rerun/FI';      % Basic name of the files ['XXXX' ('.gzip.h5')]
expnames = split(expname,'/');
disp(expnames)

setup       = 1;            % Model setup: 1=largescale, 2=lab

nstart      = 999;          % Starts from file number XXX.
nend        = nstart;       % Ends with file number XXX.
skip        = 1;            % Skip this many files. Can use nend-nstart   

% --- What variables ---
% When want to plot a single variable make (2) or (3) 0. 
sp_var(1)   = 20;           % 1 - vx ; 2 - vy ; 3 - density ; 4 - temperature ; 5 - pressure; 6 - sxx, 7 - sxy, 8 -visocisty; 9 - exx,
sp_var(2)   = 0;           % 10 - exy; 11 - G (not in all hdf5); 12 - exx_ne; 13 - exy_ne; 14 - s_ii; 15 - e_ii; 16 - e_ne_ii; 17 - Z (G not in all hdf5); 20 - composition; OWN: 21 - x-grid
sp_var(3)   = 0;           % This third never equal to 20, since colormap of last variable dictates which one use (freezeColors does not work yet)                                    
                            % v=velocity, s=stress, e=strain rate, G=shear modulus
                            
temp        = 1;            % Temperature contours

% --- How to plot it ---
printmod    = 2;        % 1: do not show on screen, but save as png (when not using GUI on servers for speed + ease! Note need equally spaced imagesc i.s.o. pcolor as has bug)
                        % 2: show and png-save figures (that can be madev_x into movies through "convert")
                        % 3: show and eps-save nice figures for papers
                        % (i.e., save as eps + rm air etc): use GUI and one file: also save data in mat for multiple plotting in different paper quality script
zoom        = 1;        % 0 or -1: no zoom = whole model;    1: zoom occurs according to limits below
axis_wrt_trench = 0;    % 0 = as simulated, 1 = trench at (0,0) 
rm_air      = 0;        % Set variable to zero in air, only works for high resolution area where grid sizes composition and var equal (so when zoom)
v_field     = 0;        % OWN ADDITION: 0 = no velocity field arrows, 1 = velocity field overlay
Slab_pull   = 0;        % OWN ADDITION: 0 = don't calculate slab pull, 1 = calculate slab pull
Drag        = 0;        % OWN ADDITION: 0 = don't calculate basal mantle drag force, 1 = calculate it    
stress_axes = 0;        % Luca dal Zilio, 2020: if 1, calculates principal stress axes. 

if Slab_pull || Drag
    timestep = 1;
    t = zeros((nend - nstart)/skip + 1, 1); % Time vector
    if Slab_pull
        SP_t = zeros((nend - nstart)/skip + 1, 1); % Slab pull vector
        A_t = zeros((nend - nstart)/skip + 1, 1); % also nice to see evolution of slab area
    end
    if Drag
        MD = zeros((nend - nstart)/skip + 1, 5, 3); % Mantle drag matrix
        PLATE_LENGTH = 700e3; % length of upper plate, hardcoded now... 
    end
end

% ------------- large-scale settings (in km) -------------
if setup==1
    
    % trench coords from comp:
    x_trench = 1050; % 840
    y_trench = 7;
    
    % Plotting limitations 
    % Papers: gxt-50 gxt+250 gyt-3 gyt+78 
    if zoom > 0
        x_beg      = x_trench-150; % 50
        x_end      = x_trench+250; % 250
        y_beg      = y_trench-3;
        y_end      = y_trench+80; % 78
    end
        
    % Isotherm location
    eq_lines = [0 100 150 350 450 1300];    %according to ref in Klingelhoefer2010. 1.9 degrees is land surface.
    
    % Sticky Air THickness from init.t3c
    SA_Th   = 20; %8 
    
    % Computational time step in s
    dt_yr = 1000; % 5
    dt = dt_yr*365.25*24*60*60;
end
%--------------------------------------------------

% Check wether should make a SubPlot with 3 variables
if sp_var(2)==0 || sp_var(3)==0
    plot3var = 0;
else
    plot3var = 1;
end

% Change colormap for Matlab R2014b
set(0,'DefaultFigureColormap',feval('copper')) % for visosity, stress
%load('/home/luuk/Documents/ETH/colormaps/ScientificColourMaps7/berlin/berlin.mat');
%set(0,'DefaultFigureColormap',flipud(berlin))
% Dock figures
if printmod==2 || printmod>3
	set(0,'DefaultFigureWindowStyle','docked'); 
else
    % Keep as large as possible for paper figure resolution
	set(0,'DefaultFigureWindowStyle','normal'); 
end

% Picture settings
save_path = sprintf('/media/luuk/My Book/MSc_paper/Figures/%s/',string(expnames(7)));
disp(save_path)
if (exist(save_path, 'dir') == 0)
    mkdir(save_path);
end
if printmod == 1
   % mkdir ./Figures/
    format = '-dpng';
    res = '-r300';
elseif printmod == 3
    format = '-deps2c';
    res = '-r600';
else
    format = '-dpng'; %'-dtiff';
    res = '-r300';
end

% 'new': new file numbering 001..010..999, 'old': 1..10..999
fntype     = 'new'; %'old';


% Catch errors
if (nstart > nend)
    error('ERROR: "nend" should be superior to "nstart"!');
end

if (printmod==1 && zoom<=0)
    %error('ERROR in plotted figure: for plotting the non-high resolution area properly you need pcolor (not imagesc) and this is only possible in GUI or with MATLAB 2015! ');
    warning('WARNING: trying to plot with pcolor instead of imagesc');
end

% ----------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------------------------



%% Loop over all desired files
for nt = nstart:skip:nend
    
    % Build filenames for loading
    if (nt<10 && strcmp(fntype, 'new') == 1)
        filename = ([expname '00' num2str(nt) '.gzip.h5']);
    elseif (nt<100 && strcmp(fntype, 'new') == 1)
        filename = ([expname '0' num2str(nt) '.gzip.h5']);
    else
        filename = ([expname num2str(nt) '.gzip.h5']);
    end
    disp(filename)

    
%     % If the file is missing open the next one
%     limit = 0;
%     while (exist(filename, 'file') == 0 && limit < 4)
%         disp(['File ', num2str(nt), ' is missing, jump to file ', num2str(nt+1)])
%         nt = nt+1;
%         if (nt<10 && strcmp(fntype, 'new') == 1)
%             filename = ([expname '00' num2str(nt) '.gzip.h5']);
%         elseif (nt<100 && strcmp(fntype, 'new') == 1)
%             filename = ([expname '0' num2str(nt) '.gzip.h5']);
%         else
%             filename = ([expname num2str(nt) '.gzip.h5']);
%         end
%         %filename = (['data/' filename]);
%         limit = limit + 1;
%     end
    
    if(exist(filename, 'file') == 0)
        %disp(['No files were found, check your data files. Exit...']);
        disp(['File ', num2str(nt), ' is missing, jump to file ', num2str(nt+1)])
        %break;
    else
    
    
    % Read the structure of the hdf5 file
    fileinfo = hdf5info(filename);
    disp(['Process file number ', num2str(nt)]);
    
    
    
    % ---------------------------------------------------------------
    % --- Load constant data once -----------------------------------
    % ---------------------------------------------------------------
    
    if nt == nstart
        x       = hdf5read(filename,'/ModelGroup/gx'); %x = x; %x = x(2:end)*1e-3;
        y       = hdf5read(filename,'/ModelGroup/gy'); y = y-SA_Th; %y = (y(2:end)*1e-3)-SA_Th;
        param   = hdf5read(filename,'/ModelGroup/Model');
        
        % Retrieves model parameters
        timesum  = param(1);
        xsize = param(2);
        ysize = param(3);
        xnumx = param(4);
        ynumy = param(5);
        
        % Build staggered nodes coordinates
        % To ensure that plot contents on right location of basic, plot them on basic node coordinates; not on staggered, central nodes coordinates
        % MATLAB already plots the coordinates you give in, in the center of its visualization cell
        x_stag  = 0.5*(x(1:end-1) + x(2:end));
        y_stag  = 0.5*(y(1:end-1) + y(2:end));
        
        % large-scale length to km
        if setup == 1
            x_stag = x_stag/1e3;
            y_stag = y_stag/1e3;
            xsize = xsize/1e3;
            ysize = ysize/1e3;
            % lab length to cm
        elseif setup == 2
            x_stag = x_stag*1e2;
            y_stag = y_stag*1e2;
            xsize = xsize*1e2;
            ysize = ysize*1e2;
        end
        
        % calculate begin and end nodes for zoom visualization with imagesc
        % since imagesc can not handle irregular grids, but pcolor can not export nice enough pictures to AI etc
        if zoom > 0
            for ix = 1:1:xnumx
                if x_stag(ix) >= x_beg
                    nx_beg = ix;
                    break;
                end
            end
            for ix = nx_beg:1:xnumx
                if (ix==xnumx)
                    nx_end = ix-1;
                    break;
                elseif x_stag(ix) >= x_end
                    nx_end = ix;
                    break;
                end
            end
            for iy = 1:1:ynumy
                if (y_stag(iy)>=y_beg)
                    ny_beg = iy;
                    break;
                end
            end
            for iy = ny_beg:1:ynumy
                if (iy==ynumy)
                    ny_end = iy-1;
                    break;
                elseif (y_stag(iy)>=y_end)
                    ny_end = iy;
                    break;
                end
            end
        else
            nx_beg = 1;
            nx_end = xnumx-1;
            ny_beg = 1;
            ny_end = ynumy-1;
        end
        
        % Correct for plotting with respect to trench
        if axis_wrt_trench == 1
            x_stag = x_stag-x_trench;
            y_stag = y_stag-y_trench-SA_Th;
            
            if setup==2
                x_topwedge = 67.1-x_trench;
                y_topwedge = 1.57-y_trench;
                
                x_trench = 0;
                y_trench = 0;
            end
        end
        
        % Make marker grid for composition
        if sp_var(1) == 20 || sp_var(2) == 20 || sp_var(3) == 20 || (rm_air==1 && setup==1)
            % Get minimum gridsize in km, since this is resolution visualization grid
            res_high = (min(y_stag(2:end)-y_stag(1:end-1)));
            
            % Get nr elements for visualization grids
            mxvislim = floor(xsize/(res_high) + 1);
            myvislim = floor(ysize/(res_high) + 1);
            
            % Composition matrix: size of dense regular markergrid
            % Whole for filling matrix correctly; zoom while plotting
            composition = ones(myvislim,mxvislim);
            
            % Calculate begin and end nodes for visualization of composition with imagesc
            % since imagesc can not handle irregular grids, but pcolor can not export nice enough pictures to AI etc
            if zoom > 0
                mnx_beg = floor(x_beg/res_high);
                mnx_end = floor(x_end/res_high);
                mny_beg = floor(y_beg/res_high);
                mny_end = floor(y_end/res_high);
            else
                mnx_beg = 1;
                mnx_end = floor(xsize/res_high);
                mny_beg = 1;
                mny_end = floor(ysize/res_high);
            end
            
            % Make marker grid
            mx = mnx_beg:1:mnx_end;
            my = mny_beg:1:mny_end;
            
            if axis_wrt_trench == 1
                mx = mx*res_high-x_trench;
                my = my*res_high-y_trench-SA_Th;
            elseif axis_wrt_trench == 0
                mx = mx*res_high;
                my = my*res_high;
            end
            
        end
        
    else
        param   = hdf5read(filename,'/ModelGroup/Model');
        timesum  = param(1);
    end
    
    
    
    % ---------------------------------------------------------------
    % --- Load variables every time step ----------------------------
    % ---------------------------------------------------------------
        
    if temp==1
        T = hdf5read(filename,'/NodeGroup/tk');
        T = reshape(T, ynumy, xnumx);         % on basic nodes
    end
    
    if plot3var == 1
        last_iv = 3;
    else
        last_iv = 1;
    end
    
    for iv = 1:1:last_iv
        vis_var = sp_var(iv);
        
        if vis_var == 1
            % make in cm/year
            var     = hdf5read(filename,'/NodeGroup/vx');
            if setup==1
                var = (var*1e2)*60*60*24*365.25;
            elseif setup==2
                var = var*100;
            end
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 2
            var     = hdf5read(filename,'/NodeGroup/vy');
            % switch around, so - is subsidence
            var = -var;
            if setup==1
                var = (var*1e2)*60*60*24*365.25;
            elseif setup==2
                var = var*100;
            end
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 3
            var     = hdf5read(filename,'/NodeGroup/ro');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 4
            var     = hdf5read(filename,'/NodeGroup/tk');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 5
            var     = hdf5read(filename,'/NodeGroup/pr');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 6
            var     = hdf5read(filename,'/NodeGroup/sxx');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 7
            var     = hdf5read(filename,'/NodeGroup/sxy');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 8
            var     = hdf5read(filename,'/NodeGroup/nu');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 9
            var     = hdf5read(filename,'/NodeGroup/exx');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 10
            var     = hdf5read(filename,'/NodeGroup/exy');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 11
            var     = hdf5read(filename,'/NodeGroup/gg');
            var = reshape(var, ynumy, xnumx);
            %exx_ne
        elseif vis_var == 12
            sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
            nd_tmp  = hdf5read(filename,'/NodeGroup/nd');
            
            sxx_tmp = reshape(sxx_tmp, ynumy, xnumx);
            nd_tmp  = reshape(nd_tmp, ynumy, xnumx);
            
            var     = sxx_tmp(1:end,1:end)./(2*nd_tmp(1:end,1:end));
            %exy_ne
        elseif vis_var == 13
            sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
            nu_tmp  = hdf5read(filename,'/NodeGroup/nu');
            
            sxy_tmp = reshape(sxy_tmp, ynumy, xnumx);
            nu_tmp  = reshape(nu_tmp, ynumy, xnumx);
            
            var     = sxy_tmp(1:end,1:end)./(2*nu_tmp(1:end,1:end));
            % sii
        elseif vis_var == 14
            sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
            sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
            
            sxx_tmp = reshape(sxx_tmp, ynumy, xnumx);
            sxy_tmp = reshape(sxy_tmp, ynumy, xnumx);
            
            sxx_tmp = sxx_tmp(1:end-1,1:end-1);     % on staggered nodes
            sxy_tmp = sxy_tmp(1:end,1:end);         % on basic nodes
            
            % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
            sxy_tmp   = 0.25*(sxy_tmp(1:end-1, 1:end-1) + sxy_tmp(2:end, 2:end) +  sxy_tmp(2:end, 1:end-1) + sxy_tmp(1:end-1, 2:end));
            
            % Calculate second invariant
            var   = sqrt(sxx_tmp.^2 + sxy_tmp.^2);
            %eii
        elseif vis_var == 15
            exx_tmp = hdf5read(filename,'/NodeGroup/exx');
            exy_tmp = hdf5read(filename,'/NodeGroup/exy');
            
            exx_tmp = reshape(exx_tmp, ynumy, xnumx);
            exy_tmp = reshape(exy_tmp, ynumy, xnumx);
            
            exx_tmp = exx_tmp(1:end-1,1:end-1);     % on staggered nodes
            exy_tmp = exy_tmp(1:end,1:end);         % on basic nodes
            
            % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
            exy_tmp   = 0.25*(exy_tmp(1:end-1, 1:end-1) + exy_tmp(2:end, 2:end) +  exy_tmp(2:end, 1:end-1) + exy_tmp(1:end-1, 2:end));
            
            % Calculate second invariant
            var   = sqrt(exx_tmp.^2 + exy_tmp.^2);
            %eii_ne
        elseif vis_var == 16
            sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
            nd_tmp  = hdf5read(filename,'/NodeGroup/nd');
            e_nexx  = sxx_tmp(1:end,1:end)./(2*nd_tmp(1:end,1:end));
            
            sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
            nu_tmp  = hdf5read(filename,'/NodeGroup/nu');
            e_nexy  = sxy_tmp(1:end,1:end)./(2*nu_tmp(1:end,1:end));
            
            e_nexx  = reshape(e_nexx, ynumy, xnumx);
            e_nexy  = reshape(e_nexy, ynumy, xnumx);
            
            e_nexx  = e_nexx(1:end-1,1:end-1);     % on staggered nodes
            e_nexy  = e_nexy(1:end,1:end);         % on basic nodes
            
            % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
            e_nexy  = 0.25*(e_nexy(1:end-1, 1:end-1) + e_nexy(2:end, 2:end) +  e_nexy(2:end, 1:end-1) + e_nexy(1:end-1, 2:end));
            
            % Calculate second invariant
            var   = sqrt(e_nexx.^2 + e_nexy.^2);
        elseif vis_var == 17
            gg     = hdf5read(filename,'/NodeGroup/gg');
            nu     = hdf5read(filename,'/NodeGroup/nu');
            
            var = (gg*dt)./(gg*dt+nu);
            
            var = reshape(var, ynumy, xnumx);           
        end
        %comp
        if vis_var == 20 || (rm_air==1 && setup==1 && iv==1)
            mtype = hdf5read(filename,'/VisMarkerGroup/Mtype');
            
            load('TMcolormap_lva','TMcolormap');
            
            % Extract composition from efficient storage
            num       = 1;
            ind       = 1;
            while num<length(mtype)
                value = mtype(num);
                
                if value==-2
                    % Compressed: the next ?? indices are given the color material
                    num_colors  =   mtype(num+1);
                    material    =   mtype(num+2);
                    ind_vec     =   ind:ind+num_colors-1;
                    
                    ind         =   ind_vec(end)+1;
                    num         =   num+3;
                elseif value==-1
                    ind_vec     =   ind;
                    material    =   0;
                    ind         =   ind_vec(end)+1;
                    num         =   num+1;
                else
                    ind_vec     =   ind;
                    material    =   value;
                    
                    ind         =   ind_vec(end)+1;
                    num         =   num+1;
                end
                
                composition(ind_vec) =   material;
            end
            
         
        end
        
        if vis_var == 21
            dx = zeros(1,xnumx-1);
            dy = zeros(1,ynumy-1);
            for i=2:xnumx
                dx(i-1) = x(i) - x(i-1);
            end
            for i=2:ynumy
                dy(i-1) = y(i) - y(i-1);
            end
            
            figure();
            subplot(1,2,1);
            plot(x);
            grid on
            subplot(1,2,2);
            plot(dx);
            figure();
            subplot(1,2,1);
            plot(y);
            grid on
            subplot(1,2,2);
            plot(dy)
        end
        % OWN ADDITION
        if v_field == 1
            vy = hdf5read(filename,'/NodeGroup/vy');
            %vy = -vy;
            vy = (vy*1e2)*60*60*24*365.25;
            vy = reshape(vy, ynumy, xnumx);

            vx = hdf5read(filename,'/NodeGroup/vx');
            vx = (vx*1e2)*60*60*24*365.25;
            vx = reshape(vx, ynumy, xnumx);      
        end 
        
        if Slab_pull == 1
            disp('Calculating slab pull')
            % Get all lithospheric mantle and Oceanic crust below LAB depth
            % get z-node of LAB
            nLAB = 1;
            [MX,MY] = meshgrid(mx, my);

            rest = abs(y_stag(nLAB) - 160);
            while (rest > 2)
                rest = abs(y_stag(nLAB) - 160);
                nLAB = nLAB + 1;
            end
%             figure(10)
%             p = pcolor(MX(:, 1:nLAB), MY(:, 1:nLAB), composition(1:end-1, 1:nLAB));
%             colormap(TMcolormap)
%             axis ij;
%             shading interp;
            % Plan: build a small loop to identify the marker that is closest
            % to each node using mx and my 
            % at that mx and my, find the composition and output it. if
            % slab, then add to slab density. Else average over
            % "surroundings" density
            var      = h5read(filename,'/NodeGroup/ro');
            rho      = reshape(var, ynumy, xnumx);
            var      = h5read(filename, '/NodeGroup/tk');
            T        = reshape(var, ynumy, xnumx);
            T_LAB    = 1617.6;
            T_sub    = T < T_LAB;
            T_subrow = sum(T_sub, 2) ; % rowwise sum of elements with sub-LAB temperatures
            [XX, YY] = meshgrid(x/1e3,y/1e3);
            rho_slab = zeros(ynumy, xnumx);
            A        = 0; % cumulative slab area
            starty   = 160e3; % depth in m below which we start calculating slab pull
            endy     = 600e3; % depth in m below which we never calculate SP
            % x and y have dimensions xnumx and ynumy
            % mx and my have dimensions mnx_end and mny_end
            for i=1:ynumy
                 if y(i) > starty % y is in meters, mx and my in km.
                     if T_subrow(i) > 100 && y(i) < endy % if less than 100 nodes have sub-LAB temperatures, the slab is not taken into account            
                         dy = abs(my - y(i)/1e3); % y difference in km
                         [mindy, minyid] = min(dy); % minimum dy and its index
                         for j=1:xnumx
                             dx = abs(mx - x(j)/1e3); % x difference in km
                             [mindx, minxid] = min(dx);% minimum dx and its index
                             dists = sqrt(mindx ^2 + mindy^2);
                             cp = composition(minyid, minxid);
                             if cp ~= 10  % if this is not asthenosphere,
                                 % this marker is slab
                                 rho_slab(i,j) = 1;                                  
                             else
                                 % this marker is no slab
                                 rho_slab(i,j) = -1;
                             end
                         end
                     else % if the slab is very thin or gone here it does not count
                        rho_slab(i:end, :) = -1; 
                        fprintf('Estimated depth of detachment %8.1f km', y(i)/1e3); 
                        break; % if this is triggered we can end the loop        
                     end
                 end
            end
            spfig = figure(11);
            set(spfig, 'Visible', 'Off');
            hax1 = axes;
            p     = pcolor(hax1,XX, YY, rho_slab);
            axis ij;
            shading interp;
            xlabel('Horizontal distance [km]');
            ylabel('Depth [km]');
            c            = colorbar;
            c.YTick      = [-1 0 1];
            c.YTickLabel = {"surroundings", "lithosphere", "slab"};
            set(gca, 'DataAspectRatio', [1 1 1]);
            title(['Time = ', num2str(timesum/1e6),' Myrs for model ',char(expnames(7))]);
            print(format, res,[save_path, char(expnames(7)),'_slab_id', num2str(nt)]);
            close(spfig);
            % Actual slab pull calculation
            % we know now where the slab is.
            SP = 0; % Total slab pull in N/m
            ystep = y(2:end) - y(1:end-1) ; % dy vector
            for i=1:ynumy-1
                r   = rho_slab(i,:); % get row of slab identification matrix
                idx = find(r > 0);   % find slab indices
                if any(idx > 0)
                    % calculate incremental area and slab pull
                    dA = ystep(i) * (x(idx(end)) - x(idx(1))); % get dA
                    rho_m = rho(i, r<0) ; 
                    rho_s = rho(i, r>0) ;
                    SP = SP + 9.81 * dA * (mean(rho_m) - mean(rho_s));
                    A = A + dA;
                end
            end
            t(timestep) = timesum;
            SP_t(timestep) = SP;
            A_t(timestep) = A;
        end
        
         if Drag
            disp('Calculating drag force')
         
            nLAB    = 1;
            [MX,MY] = meshgrid(mx, my);
            rest    = abs(y_stag(nLAB) - 160);
            while (rest > 2)
                rest = abs(y_stag(nLAB) - 160);
                nLAB = nLAB + 1;
            end
%           % identify x_start by looking at the non-asthenosphere node
%           from the right
            var     = h5read(filename, '/NodeGroup/vx');
            vx      = reshape(var, ynumy, xnumx);
            dVxdz   = vx(2:end, 2:end) ./ (y(2:end) - y(1:end-1)); % vertical gradient of horizontal velocity
            var     = h5read(filename, '/NodeGroup/nu');
            eta_eff = reshape(var, ynumy, xnumx);
            depths  = nLAB-2:1:nLAB+2;
           
            % get x_start and x_end: get first non-asthenosphere node 
            nUP_start = 0; % node at which upper plate starts
            % very similar to SP calculation... 
            for i=1:ynumy
                 if i > nLAB && i <= nLAB + 2 % y is in meters, mx and my in km.
                     dy = abs(my - y(i)/1e3); % y difference in km
                     [mindy, minyid] = min(dy); % minimum dy and its index
                     for j=1:xnumx
                         dx = abs(mx - x(j)/1e3); % x difference in km
                         [mindx, minxid] = min(dx);% minimum dx and its index
                         dists = sqrt(mindx ^2 + mindy^2);
                         cp = composition(minyid, minxid);
                         if cp ~= 10  
                             % this marker is slab
                             nUP_start = j;  
                             break;
                         end
                     end
                     if nUP_start > 0
                         break;
                     end
                 end
            end
            if nUP_start > 0 % if any non-asthenoshperic node was found
                nUP_end = find(x > (x(nUP_start) + PLATE_LENGTH),1, 'first');
                %disp(["using ", x(nUP_end)-x(nUP_start), " m for the integration"])
                if nUP_end > nUP_start
                    % get all components for integration of shear stress,
                    % which is roughly the same as integration of the
                    % product of dVxdz and the effective viscosity
                    dx        = x(2:end) - x(1:end-1);
                    dVxdz     = dVxdz(depths, nUP_start:nUP_end);
                    eta_eff   = eta_eff(depths, nUP_start:nUP_end);
                    tractions = dVxdz .* eta_eff;
                    R         = repmat((dx(nUP_start:nUP_end))', 5, 1);
                    Dragforce = tractions .* R;            
            
                    % store results in matrix before next time step
                    t(timestep)        = timesum;
                    MD(timestep, :, 1) = mean(-Dragforce, 2);
                    MD(timestep, :, 2) = max(-Dragforce,[], 2);
                    MD(timestep, :, 3) = min(-Dragforce,[], 2);
                end
            else
                disp('no non-asthenosphere nodes found')
            end
        end
        % END OWN ADDITION

        %% Remove velocities in air for paper plotting to avoid distraction
        if rm_air==1
            if setup==1
                if zoom>0
                    for i=1:1:length(mx)
                        for j=1:1:length(my)
                            % adjust composition and var arrays such that cover same area
                            if composition(mny_beg+j,mnx_beg+i)==0
                                var(ny_beg+j,nx_beg+i) = 1e-10;
                            end
                        end
                    end
                else
                    disp('NOTE: Can not remove air for no zoom (need all within HR grid): will never do for papers anyway!');
                end
                
            elseif setup == 2
                % Calculate y-coordinate of top surface wedge for each x (with a bit cut off)
                slope = (y_topwedge-y_trench)/(x_topwedge-x_trench);
                begin = y_trench - slope*x_trench;
                
                % Set following part to near 0 (not 0 for if strain rates)
                for i=1:1:xnumx-1
                    for j=1:1:ynumy-1
                        lwsurface = slope*x_stag(i)+begin;
                        if (x_stag(i) < x_trench) && (y_stag(j) < y_trench)
                            var(j,i) = 1e-10;
                        elseif (y_stag(j) < y_topwedge)
                            var(j,i) = 1e-10;
                        elseif (x_stag(i) > x_trench) && (x_stag(i) < x_topwedge) && (y_stag(j) < lwsurface)
                            var(j,i) = 1e-10;
                        end
                    end
                end
                % note: neglected top of wall; made it zero as well, is not moving anyway
            end
        end
        
        % Fill corresponding variable for subplotting 3 variables 
        if plot3var == 1 
            if iv==1 && sp_var(iv)~=20
                var1 = var;
            elseif iv==2 && sp_var(iv)~=20
                var2 = var;
            elseif iv==3 && sp_var(iv)~=20
                var3=var;
            end
        end
    end
        
    
    % ---------------------------------------------------------------
    % --- Plot figure -----------------------------------------------
    % ---------------------------------------------------------------
    
    if printmod == 1
        fig=figure(1)  ;  
        set(fig, 'Visible', 'Off');
    elseif printmod == 3
        figure('units','normalized','outerposition',[0 0 1 1]);
    else
        figure(1);
        clf;
    end
    
    for iv = 1:1:last_iv
        
        % Load appropriate variables for this iv
        vis_var = sp_var(iv);
        
        if plot3var == 1
            subplot(3,1,iv);
            
            if iv == 1 && sp_var(iv)~=20
                var = var1;
            elseif iv == 2 && sp_var(iv)~=20
                var = var2;
            elseif iv == 3 && sp_var(iv)~=20
                var = var3;
            end
        end
        
        % Log variables
        if vis_var==8 || vis_var==14 || vis_var==15 || vis_var==16
            colormap('default');
            if printmod==2 || printmod>3
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), log10(var(ny_beg:ny_end,nx_beg:nx_end)));
            elseif printmod==1 || printmod==3
               % imagesc(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), log10(var(ny_beg:ny_end,nx_beg:nx_end)));		
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), log10(var(ny_beg:ny_end,nx_beg:nx_end)));
               
            end
            %freezeColors;
            cb=colorbar;
            % Composition
        elseif vis_var == 20
            if plot3var == 0
                colormap(TMcolormap);
            end
            if axis_wrt_trench == 1
                imagesc(mx, my, composition(mny_beg:1:mny_end,mnx_beg:1:mnx_end));
            elseif axis_wrt_trench == 0
                imagesc(mx, my, composition(mny_beg:1:mny_end,mnx_beg:1:mnx_end));
            end
      
            % daspect([1 2 1]);
            %freezeColors;
            cb=colorbar;
            % All other variables
        else
            colormap('default');
            if printmod > 1
                pcolor(x_stag(nx_beg:5:nx_end), y_stag(ny_beg:5:ny_end), var(ny_beg:5:ny_end,nx_beg:5:nx_end));
                colormap(redblue)
            else
                %imagesc(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), var(ny_beg:ny_end,nx_beg:nx_end));
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), var(ny_beg:ny_end,nx_beg:nx_end));
            end
            %freezeColors;
            cb=colorbar;
        end
        % OWN ADDITION: wanted to overlay composition with velocity field,
        % but never found the time to get it to work nicely with the
        % irregular grid
        if v_field == 1
            step = 10;
            vskip = 40; % vertical nodes to skip (to prevent air velocities)
            hold on;
            [X,Y]       = meshgrid(x_stag(nx_beg:step:nx_end), y_stag(ny_beg + vskip:step:ny_end));
            q           = quiver(X,Y,vx(ny_beg + vskip:step:ny_end,nx_beg:step:nx_end),vy(ny_beg + vskip:step:ny_end,nx_beg:step:nx_end));
            %q.AutoScale = 'off';
            q.Color     = 'black';
%             hold on;
%             quiver(2750, 780, 20, 0);
%             text(2850, 760, '1 cm/yr')
        end
            % END OWN ADDITION
        
                    
  %============================================
        % Stress axes (Luca dal Zilio, 2020)
        if stress_axes == 1
            %============================================
            var_sxx = hdf5read(filename,'/NodeGroup/sxx');
            sxx = reshape(var_sxx, ynumy, xnumx);
            syy = -sxx;
            var_sxy = hdf5read(filename,'/NodeGroup/sxy');
            sxy = reshape(var_sxy, ynumy, xnumx);
            
            for m=2:ynumy-1
                for n=2:xnumx-1
                    SXX(m,n)=(sxx(m,n)+sxx(m+1,n)+sxx(m,n+1)+sxx(m+1,n+1))/4;
                    SYY(m,n)=(syy(m,n)+syy(m+1,n)+syy(m,n+1)+syy(m+1,n+1))/4;
                    sxy(m,n)=sxy(m,n);
                    sii(m,n)= sqrt(SXX(m,n)^2 + sxy(m,n)^2);
                end
            end
 
            for m=40:13:160 % MODIFY: Y 
                for n=260:26:1040-1 % MODIFY: X
                    Stress=[SXX(m,n) sxy(m,n);sxy(m,n) SYY(m,n)];
                    [a,b]=eig(Stress);
                    if (b(1,1)>b(2,2))
                        stress_max=b(1,1);
                        column_max=1;
                        stress_min=b(2,2);
                        column_min=2;
                    else
                        stress_max=b(2,2);
                        column_max=2;
                        stress_min=b(1,1);
                        column_min=1;
                    end
                    %===============================================
                    stress_max=stress_max/sii(m,n);
                    stress_min=abs(stress_min)/sii(m,n);
                    %Plot max and min axes
                    hold all
                    %===============================================
                    %Maximum extension: black
                    %vector length and scale it
                    l=-stress_max:0.01:stress_max;
                    l=l*2;
                    plot(x_stag(n)+l*a(1,column_max),y_stag(m)+l*a(2,column_max),'k','LineWidth',2)
                    %===============================================
                    %Maximum compression: red
                    %vector length and scale it
                    l=-stress_min:0.01:stress_min;
                    l=l*2;
                    plot(x_stag(n)+l*a(1,column_min),y_stag(m)+l*a(2,column_min),'color','r','LineWidth',2)
                end
            end
        end
        
  %====================================================
            
        % Isotherms
        if temp == 1
            hold on;
            [C2,h2] = contour(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), T(ny_beg:ny_end,nx_beg:nx_end)-273.15,eq_lines,'w','Linewidth',0.5);
            axis ij image;
            hold off;
        end        
        
        % ---------------------------------------------------------------
        % --- Set visual properties -------------------------------------
        % ---------------------------------------------------------------
        
        % --- large-scale settings ---
        if setup==1
            
            if vis_var == 1
                ylabel(cb,'Horizontal velocity [cm/yr]')
                vx_lim = 1.5e-9;
                caxis ([-1 1]); % [-1 6]
            elseif vis_var == 2
                ylabel(cb,'Vertical velocity [cm/yr]')
                vy_lim = 1e-9;
                caxis ([-1 1]); %[-3 3]
            elseif vis_var == 3
                ylabel(cb,'Density [kg/m^3]');
                caxis ([2500 3350]);
            elseif vis_var == 4
                ylabel(cb,'Temperature [dK]');
                caxis ([750 850]);
                %caxis ([1500 1700]); % [273 1700]
            elseif vis_var == 5
                ylabel(cb,'Pressure [Pa]');
                %caxis ([1.3e9 1.4e9]);
            elseif vis_var == 6
                ylabel(cb,'Dev. normal stress [Pa]');
                caxis ([-50e6 50e6]);
            elseif vis_var == 7
                ylabel(cb,'Dev. shear stress [Pa]');
                caxis ([-30e6 30e6]);
            elseif vis_var == 8
                caxis ([17 25]);
                ylabel(cb,'log10(Viscosity) [Pa s]');
            elseif vis_var == 9
                ylabel(cb,'\epsilon_{xx} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 10
                ylabel(cb,'\epsilon_{xy} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 11
                ylabel(cb,'shear modulus [Pa]');
                caxis ([1e10 8e10]);
            elseif vis_var == 12
                ylabel(cb,'\epsilon_{vp, xx} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 13
                ylabel(cb,'\epsilon_{vp, xy} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 14
                ylabel(cb,'log10(stress II) [Pa]');
                caxis(gca,[6 8]);
            elseif vis_var == 15
                ylabel(cb,'log10(strain rate II) [s^{-1}]');
                caxis(gca,[-15 -12]);
            elseif vis_var == 16
                ylabel(cb,'log10(ne strain rate II) [s^{-1}]');
                caxis(gca,[-15 -12]);
            elseif vis_var == 17
                ylabel(cb,'Visco-elasticity factor Z [-]');
                caxis(gca,[0 1]);
            elseif vis_var == 20
                ylabel(cb,'rock composition number');
                if plot3var == 1
                    caxis ([0 19]);
                else
                    caxis ([0 19]);
                end
            end
            
            % --- lab settings ---
        elseif setup==2
            
            if vis_var == 1
                ylabel(cb,'Horizontal velocity [cm/s]')
                caxis ([-0.05 0.005]);
            elseif vis_var == 2
                ylabel(cb,'Vertical velocity [cm/s]')
                caxis ([-0.5e-3 1e-3]);
            elseif vis_var == 3
                ylabel(cb,'Density [kg/m^3]');
                caxis(gca,[900 1200])
            elseif vis_var == 4
                ylabel(cb,'Temperature [dK]');
                %caxis ([400 700]);
            elseif vis_var == 5
                ylabel(cb,'Pressure [Pa]');
                %caxis ([1.3e9 1.4e9]);
            elseif vis_var == 6
                ylabel(cb,'Dev. normal stress [Pa]');
                %caxis ([-50e6 50e6]);
            elseif vis_var == 7
                ylabel(cb,'Dev. shear stress [Pa]');
                %caxis ([-30e6 30e6]);
            elseif vis_var == 8
                caxis ([0 8]); %2.7 3.2
                ylabel(cb,'log10(Viscosity) [Pa s]');
            elseif vis_var == 9
                ylabel(cb,'\epsilon_{xx} [s^{-1}]');
                %caxis ([-2e-14 2e-14]);
            elseif vis_var == 10
                ylabel(cb,'\epsilon_{xy} [s^{-1}]');
                %caxis ([-2e-14 2e-14]);
            elseif vis_var == 11
                ylabel(cb,'shear modulus [Pa]');
                caxis ([1e10 1e11]);
            elseif vis_var == 12
                ylabel(cb,'\epsilon_{vp, xx} [s^{-1}]');
                %caxis ([-2e-14 2e-14]);
            elseif vis_var == 13
                ylabel(cb,'\epsilon_{vp, xy} [s^{-1}]');
                %caxis ([-2e-14 2e-14]);
            elseif vis_var == 14
                ylabel(cb,'log10(stress II) [Pa]');
                caxis ([1.6 2.2]);
            elseif vis_var == 15
                ylabel(cb,'log10(strain rate II) [s^{-1}]');
                %caxis(gca,[-15 -12]);
            elseif vis_var == 16
                ylabel(cb,'log10(ne strain rate II) [s^{-1}]');
                %caxis(gca,[-15 -12]);
            elseif vis_var == 17
                ylabel(cb,'Visco-elasticity factor Z [-]');
                caxis(gca,[0 1]);
            elseif vis_var == 20
                ylabel(cb,'rock composition number');
                caxis ([0 19]);
            end
        end
        
        if vis_var==20
            shading flat
        else
            shading interp
        end
        axis ij image
        if (plot3var == 0 || iv == 3)
            xlabel('X');
        end
        ylabel('Z');
        box on;
        set(gca, 'Layer','top');
        %set(gca,'xTick',750:100:1150);
        %set(gca,'yTick',0:50:200);
        
        % Add title
	% subtract path from expname
	final_expname = char(expnames(7));
        if (iv == 1)
            if setup == 1
                title(['Time = ', num2str(timesum/1e6),' Myrs for model ',final_expname]);
            elseif setup == 2
                title(['Time = ', num2str(timesum*cor_sd),' s for model ',final_expname]);
            end
        end
    end

    
    
    % ---------------------------------------------------------------
    % --- Save figures ----------------------------------------------
    % ---------------------------------------------------------------
    
    if printmod >= 1
        % Build up filename
        if vis_var == 1
            fig_var = 'vx';
        elseif vis_var == 2
            fig_var = 'vy';
        elseif vis_var == 3
            fig_var = 'rho';
        elseif vis_var == 4
            fig_var = 'tk';
        elseif vis_var == 5
            fig_var = 'pr';
        elseif vis_var == 6
            fig_var = 'sxx';
        elseif vis_var == 7
            fig_var = 'sxy';
        elseif vis_var == 8
            fig_var = 'nu';
        elseif vis_var == 9
            fig_var = 'exx';
        elseif vis_var == 10
            fig_var = 'exy';
        elseif vis_var == 11
            fig_var = 'gg';
        elseif vis_var == 12
            fig_var = 'exx_ne';
        elseif vis_var == 13
            fig_var = 'exy_ne';
        elseif vis_var == 14
            fig_var = 'sii';
        elseif vis_var == 15
            fig_var = 'eii';
        elseif vis_var == 16
            fig_var = 'eii_ne';
        elseif vis_var == 17
            fig_var = 'Z';
        elseif vis_var == 20
            fig_var = 'comp';
            if v_field == 1
                fig_var = append(fig_var,'_v_');
            end
%             if Slab_pull
%             fig_var = 's_id';
%             end
        elseif vis_var == 21
            fig_var = 'xgrid';
        
        end
        
        
        if zoom > 0 && plot3var == 0
            ext = '_z_';
            set(gca, 'BoxStyle', 'full', 'Box', 'on','DataAspectRatio',[1 1 1]);
        elseif zoom > 0 && plot3var == 1
            ext = '_zsp_';
        elseif zoom <= 0 && plot3var == 0
            ext = '_';
        elseif zoom <= 0 && plot3var == 1
            ext = '_sp_';
        end
        
        if nt<10
            filename_fig    =  [save_path,final_expname,ext,fig_var,'00', num2str(nt)];
        elseif nt<100 &&  nt>=10
            filename_fig    =  [save_path, final_expname,ext,fig_var,'0', num2str(nt)];
        else
            filename_fig    =  [save_path, final_expname,ext,fig_var, num2str(nt)];
        end
        % Note adds last variable to filename in case of subplot
                
        % Can uncomment this print comment if testing many figures
        print(format,res,filename_fig);
        % For paper quality figures:
        % Note imagesc and deps2c allow you to nicely edit figures in Illustrator (prior tp CS6)
        
        % Data save so that can very nice plots in a different script
        % if desired (if printmod=3 quality or complexity is not enough)
        if printmod==3 && nt==nstart
            if zoom <= 0
                file_save = ['hq_fig_data_',fig_var,'_',char(expname),num2str(nt),'.mat'];
            else
                file_save = ['hq_fig_data_z_',fig_var,'_',char(expname),num2str(nt),'.mat'];
            end
            if vis_var==20
                if setup == 1
                    save(file_save,'mx','my','composition','T','x','y','x_stag','y_stag','res_high','y_trench','x_trench','SA_Th','nx_beg','ny_beg','nx_end','ny_end','mnx_beg','mny_beg','mnx_end','mny_end','eq_lines');
                elseif setup == 2
                    save(file_save,'mx','my','composition','x','y','x_stag','y_stag','res_high','y_trench','x_trench','nx_beg','ny_beg','nx_end','ny_end','mnx_beg','mny_beg','mnx_end','mny_end');
                end
            else
                if setup == 1
                    save(file_save,'var','T','x','y','x_stag','y_stag','y_trench','x_trench','SA_Th','nx_beg','ny_beg','nx_end','ny_end','eq_lines');
                elseif setup == 2
                    save(file_save,'var','x','y','x_stag','y_stag','y_trench','x_trench','nx_beg','ny_beg','nx_end','ny_end');
                end
            end
        end   
        disp(['Image ', num2str(nt), ' printed and saved.'])
        
        if Slab_pull
            fprintf('Time %f,Slabpull of %f N, area %f m^2\n', timesum/1e6, SP, A)
        end
        if Slab_pull || Drag
            timestep = timestep + 1;
        end
        
%         if Slab_pull
%             D = [t SP_t A_t];
%             disp('writing Slabpull data to file')
%             fn = sprintf('%s_SPdata.txt', final_expname);
%             save(fn, 'D', '-ascii', '-double', '-tabs');
%         end
%         if Drag
%             disp('writing Suction data to files')
%             fn_mean = sprintf('%s_MDdata_mean.txt', final_expname);
%             fn_max = sprintf('%s_MDdata_max.txt', final_expname);
%             fn_min = sprintf('%s_MDdata_min.txt', final_expname);
%             MDmean = [t MD(:,:,1)];
%             MDmax = [t MD(:,:,2)];
%             MDmin = [t MD(:,:,3)];
%             save(fn_mean, 'MDmean', '-ascii', '-double', '-tabs');
%             save(fn_max, 'MDmax', '-ascii', '-double', '-tabs');
%             save(fn_min, 'MDmin', '-ascii', '-double', '-tabs');
% 
%         end
    end
    end
    
end
% OWN ADDITION: Slab pull and mantle drag data saving to files
if Slab_pull
    D = [t SP_t A_t];
    disp('writing Slabpull data to file')
    fn = sprintf('%s_SPdata.txt', final_expname);
    save(fn, 'D', '-ascii', '-double', '-tabs');
end
if Drag
    disp('writing Suction data to files')
    fn_mean = sprintf('%s_MDdata_mean.txt', final_expname);
    fn_max = sprintf('%s_MDdata_max.txt', final_expname);
    fn_min = sprintf('%s_MDdata_min.txt', final_expname);
    MDmean = [t MD(:,:,1)];
    MDmax = [t MD(:,:,2)];
    MDmin = [t MD(:,:,3)];
    save(fn_mean, 'MDmean', '-ascii', '-double', '-tabs');
    save(fn_max, 'MDmax', '-ascii', '-double', '-tabs');
    save(fn_min, 'MDmin', '-ascii', '-double', '-tabs');
    
end
if nt~=nend
    close(gcf);
end   



