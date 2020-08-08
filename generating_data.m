%% This generates the data for the different numerical phantoms with different noise levels.


signal_to_noise_ratio1=60;

%clear all;
% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
% create the time array
time.dt = 5e-8;         % sampling time in sec 5e-8
time.length = 500;      % number of points in time

object_sim.Nx = 501;  % number of grid points in the x (row) direction
object_sim.Ny = 501;  % number of grid points in the y (column) direction
object_sim.x = 50.1e-3;              % total grid size [m]
object_sim.y = 50.1e-3;              % total grid size [m]

dx = object_sim.x/object_sim.Nx;              % grid point spacing in the x direction [m]
dy = object_sim.y/object_sim.Ny;              % grid point spacing in the y direction [m]
kgrid = makeGrid(object_sim.Nx, dx, object_sim.Ny, dy);

% create a second computation grid for the reconstruction to avoid the
% inverse crime
object_rec.Nx = 1001;  % number of grid points in the x (row) direction
object_rec.Ny = 1001;  % number of grid points in the y (column) direction
object_rec.x = 50.1e-3;              % total grid size [m]
object_rec.y = 50.1e-3;              % total grid size [m]

dx_rec = object_rec.x/object_rec.Nx;              % grid point spacing in the x direction [m]
dy_rec = object_rec.y/object_rec.Ny;              % grid point spacing in the y direction [m]
recon_grid = makeGrid(object_rec.Nx, dx_rec, object_rec.Ny, dy_rec);

% define a centered Cartesian circular sensor
clear sensor;
sensor_radius = 22e-3;     % [m]
sensor_angle = 2*pi;      % [rad]
sensor_pos = [0, 0];        % [m]
num_sensor_points = 100;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points, sensor_pos, sensor_angle);
sensor.mask = cart_sensor_mask;
sensor.frequency_response = [2.25e6 70];



M = 100;
N = 100;
indxi = ceil(object_sim.Nx/2) - M:ceil(object_sim.Nx/2) +M;
indyi = ceil(object_sim.Ny/2) -N:ceil(object_sim.Ny/2) +N;

Nxi = length(indxi);
Nyi = length(indyi);
tl = time.length;
sml = length(sensor.mask);
ANx = tl*sml;
ANy = Nxi*Nyi;
Nx = tl*sml;
object=object_sim;
xx =-100;
yy =-100;
c_x = ceil(object.Nx/2); c_y = ceil(object.Ny/2);


signal_to_noise_ratio = signal_to_noise_ratio1;    % [dB]

%% phantom blood vessel

A=imread('my_pm7.png');
scale=0.5;
B = imresize(A,[420 440]);
imshow(B,[]);

% bv = imread('my_pm6.png');
bv=B;
BV2=im2bw(bv,0.9);
in1=find(BV2==1);
in0=find(BV2==0);
BV2(in1)=1;
BV2(in0)=0;

BV2 = (BV2(5:405,20:420));

Mr = 200;
Nr = 200;
indxr = ceil(object_rec.Nx/2) - Mr:ceil(object_rec.Nx/2) +Mr;
indyr = ceil(object_rec.Ny/2) -Nr:ceil(object_rec.Ny/2) +Nr;
object_rec.p0 = zeros(object_rec.Nx, object_rec.Ny);
object_rec.p0(indxr,indyr) = BV2(:,:);
sensor.frequency_response = [2.25e6 70];
sd2 = forward(object_rec, time, medium, sensor);
% add noise to the recorded sensor data
sdn2 = addNoise(sd2, signal_to_noise_ratio, 'peak');
sdn2_v_bv=reshape(sdn2,50000,1);


for i = 1:401
    it1 = BV2(i,1:2:end);
    BVf(i,:) = it1(:);
end
BV2_bv = BVf(1:2:end,:);
% figure,imshow(BV2_bv,[]);colorbar;

%% derenzo phantom

bv = imread('derenzo.png');
bv = double(rgb2gray(bv));
ind = find(bv(:)<120);
bv(ind) = 0;
ind = find(bv(:)>120);
bv(ind)=1;
BV1 = imresize(double(bv),[430 430]);
ind = find(BV1(:)<0.9);
BV1(ind) = 0;
ind = find(BV1(:)>0.9);
BV1(ind) = 1;

BV2 = BV1(10:410,20:420);

Mr = 200;
Nr = 200;
indxr = ceil(object_rec.Nx/2) - Mr:ceil(object_rec.Nx/2) +Mr;
indyr = ceil(object_rec.Ny/2) -Nr:ceil(object_rec.Ny/2) +Nr;
object_rec.p0 = zeros(object_rec.Nx, object_rec.Ny);
object_rec.p0(indxr,indyr) = BV2(:,:);
sensor.frequency_response = [2.25e6 70];
sd2 = forward(object_rec, time, medium, sensor);
% add noise to the recorded sensor data
sdn2 = addNoise(sd2, signal_to_noise_ratio, 'peak');
sdn2_v_der=reshape(sdn2,50000,1);


for i = 1:401
    it1 = BV2(i,1:2:end);
    BVf(i,:) = it1(:);
end
BV2_der = BVf(1:2:end,:);

% imshow(reshape(A_b'*sdn2_v,201,201),[]);