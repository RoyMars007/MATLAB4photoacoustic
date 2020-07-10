clc
clear all
%==================================================%
% Press_distribute: the initial distribute of A-line collection          %
% dy: element spacing; dt(x): flight time (displacement)              %
% c: wave velocity; img_scale: 1e3-->[mm];1e6-->[um]        %
%delay: To remove the trigger signal, it is OK if you do not know.%
% This is for PA linear-array reconstruction in 2D (3D)                %
% The kWave's built-in class has been employed in this method  %
% Roy Ma July, 10th, 2020                                                           %
%==================================================%
% start timer
tic;
%====================defination=========================%
num_req_inputs = 4;
interp_method = '*linear';%'linear' (default) | 'nearest' | 'cubic' | 'spline' | 'makima'
img_scale = 1e3; % enlarged image with unit [mm]
dy = (0.3)*1e-3; % element spacing 0.3 [mm]
dt = 1/18*1e-6;  % 18 [MHz]
c = 1580; % velocity of PA wave
delay=101; % tuncation trigger
tunc_length = 1896;

% read data 2D, 
% For each B-scan, we collect the signal with the cycle 128 columns per time and collect total 20 times.
temp = dlmread('D:\YUAN\2\Slot11.txt');
L = size(temp,1);
W = size(temp,2);
temp_sum = zeros(L/20,W);
for idx=1:L/128
    temp_sum = temp_sum + temp((1+(idx-1)*128):idx*128,:);
end
press_distribute = temp_sum./20;
press_distribute = press_distribute(:,delay:tunc_length)';%yt-->ty


%======================calculation=======================%
% reversed order of p
press_distribute = [flip(press_distribute, 1); press_distribute(2:end, :)];
[Nt, Ny] = size(press_distribute);
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
press_distribute = sf .* fftshift(fftn(ifftshift(press_distribute)));
% condition
press_distribute(abs(w) < abs(c * comsol_grid.ky)) = 0;
press_distribute = interp2(comsol_grid.ky, w, press_distribute, comsol_grid.ky, w_new, interp_method);
press_distribute(isnan(press_distribute)) = 0;
% compute the inverse FFT of p(kx, ky) to yield p(x, y)
press_distribute = real(fftshift(ifftn(ifftshift(press_distribute))));
% remove the left part of the mirrored data which corresponds to the negative part of the mirrored time data
press_distribute = press_distribute( ((Nt + 1) / 2):Nt, :);

% correct the scaling the forward FFT, dy = dt*c,  only half the plane collects (symmetrical£¬the limited view problem)
press_distribute = 2 * 2 * press_distribute ./ c;

% end timer
toc
%================ the specified imagesc ================%
max_pd = max(press_distribute(:));
x_axis = [0, (Nt / 2) * dt * c];
y_axis = [0, Ny * dy];
figure,imagesc('XData',y_axis.*img_scale,'YData',x_axis.*img_scale,'CData',press_distribute, [0,max_pd]), colormap(hot);
axis image
if img_scale ==1e3
    xlabel(['Sensor Position [mm]']);
    ylabel(['Depth [mm]']);
elseif img_scale ==1e6
    xlabel(['Sensor Position [um]']);
    ylabel(['Depth [um]']);
end
    