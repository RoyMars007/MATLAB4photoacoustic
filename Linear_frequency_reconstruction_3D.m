
clc
clear all
%==================================================%
% Press_distribute: the initial distribute of A-line collection          %
% dy: element spacing; dt(x): flight time (displacement);              %
% dz=B-scan spacing                                                                     %              
% c: wave velocity; img_scale: 1e3-->[mm];1e6-->[um]             %
%delay: To remove the trigger signal, it is OK if you do not know. %
% This is for PA linear-array reconstruction in 2D (3D)                 %
% The kWave's built-in class has been employed in this method   %
% Roy Ma July, 10th, 2020                                                             %
%================================================== %
% start timer
tic;
%====================defination=======================%
num_req_inputs = 4;
interp_method = '*linear';%'linear' (default) | 'nearest' | 'cubic' | 'spline' | 'makima'
img_scale = 1e3; % enlarged image with unit [mm]
dy = (0.3)*1e-3; % element spacing 0.3 [mm]
dt = 1/18*1e-6;  % 18 [MHz]
dz = (0.3)*1e-3;
c = 1580; 
delay=101; 
tunc_length = 1896;
Nz = 4;

%=================read data from txt===================%
% For each B-scan, we collect the signal with the cycle 128 columns
% per time and collect total 20 times.
for Bscan_count = 1:Nz
    path = ['D:\YUAN\2\Slot1' num2str(Bscan_count) '.txt'];
    temp = dlmread(path);
    L = size(temp,1);
    W = size(temp,2);
    temp_sum = zeros(L/20,W);
    for idx=1:L/128
        temp_sum = temp_sum + temp((1+(idx-1)*128):idx*128,:);
    end
    pressure_distribute = temp_sum./20;
    pressure_distribute = pressure_distribute(:,delay:tunc_length)';%yt-->ty


    %====================calculation======+==============%
    % reversed order of p
    pressure_distribute = [flip(pressure_distribute, 1); pressure_distribute(2:end, :)];
    [Nt, Ny] = size(pressure_distribute);
    % KwaveGrid is utilized to construct calculation grid
    comsol_grid = comsolGrid(Nt, dt * c, Ny, dy);

    % dt*c=dx means the x/t for PA wave axis direction
    w = c .* comsol_grid.kx;% kx = kt, that means w = c*dt = x direction (axis direction's displacement)£¬

    % remap the computational grid for kx onto w, w/c = (kx^2 + ky^2)^1/2 = k. This is used for the interpolation from p(w, ky)
    % to p(kx, ky). Only real w is taken to force kx (and thus x) to be symmetrical about 0 after the interpolation. 
    w_new = c .* comsol_grid.k; % s/t, displacement
    % calculate the scaling factor using the value of kx, where kx = sqrt( (w/c).^2 - kgrid.ky.^2 ) and then manually
    sf = c.^2 .* sqrt( (w ./ c).^2 - comsol_grid.ky.^2) ./ (2 .* w);
    sf(w == 0 & comsol_grid.ky == 0) = c ./ 2;
    % compute the FFT of the input data p(t, y) to yield p(w, ky) and scale
    pressure_distribute = sf .* fftshift(fftn(ifftshift(pressure_distribute)));
    % condition
    pressure_distribute(abs(w) < abs(c * comsol_grid.ky)) = 0;
    pressure_distribute = interp2(comsol_grid.ky, w, pressure_distribute, comsol_grid.ky, w_new, interp_method);
    pressure_distribute(isnan(pressure_distribute)) = 0;
    % compute the inverse FFT of p(kx, ky) to yield p(x, y)
    pressure_distribute = real(fftshift(ifftn(ifftshift(pressure_distribute))));
    % remove the left part of the mirrored data which corresponds to the negative part of the mirrored time data
    pressure_distribute = pressure_distribute( ((Nt + 1) / 2):Nt, :);

    % correct the scaling the forward FFT, dy = dt*c,  only half the plane collects (symmetrical£¬the limited view problem)
    pressure_distribute = 2 * 2 * pressure_distribute ./ c;
    % end timer in loop
    toc

    %================ the specified imagesc ================%
    max_pd = max(pressure_distribute(:));
    x_axis = [0, (Nt / 2) * dt * c];
    y_axis = [0, Ny * dy];
    figure(1),imagesc('XData',y_axis.*img_scale,'YData',x_axis.*img_scale,'CData',pressure_distribute, [0,max_pd]), colormap(hot);
    axis image
    if img_scale ==1e3
        xlabel(['Sensor Position [mm]']);
        ylabel(['Depth [mm]']);
    elseif img_scale ==1e6
        xlabel(['Sensor Position [um]']);
        ylabel(['Depth [um]']);
    end
    title(['B-scan:' num2str(Bscan_count)]);
    pressure_3D(Bscan_count, :, :) = pressure_distribute;
end

%========================MIP==========================%
MIP = squeeze(max(pressure_3D,[],2));
figure(2),
z_axis = [0,Nz *dz];
max_mip = max(MIP(:));
imagesc('XData',y_axis.*img_scale, 'YData',z_axis.*img_scale, 'CData', MIP, [0,max_mip]), colormap(hot);
axis image

% end timer, finally...
toc
