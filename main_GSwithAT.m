%% Code for performing GS phase retrieval aided wth automatic affine 
% transform in digital in-line holographic microscopy
% 
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
%
% Last modified: 05.01.2024

clear; close all; clc

% Load holograms
for tt = 1:6
    H(:,:,tt) = double(imread(['./Data/H',num2str(tt),'.png']));
end

% System parameters
Z = [11.4, 11.87, 12.4, 12.93, 13.45, 14]*1e3; % propagation distances [um]
lambda = 0.405; % wavelength [um]
dx = 2.4; % effective pixel size [um]
n0 = 1; % refractive index of medium (air)
%% Display holograms (2 ROIs)
% close all
x1 = 1001:1600; y1 = 851:1250;
x2 = 4001:4600; y2 = 1201:1600;
for tt = 1:size(H,3)
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    tiledlayout('flow')
    nexttile; imagesc(H(y1,x1,tt),[0,255]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 1']); axis off
    nexttile; imagesc(H(y2,x2,tt),[0,255]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 2']); axis off
    pause(0.1)
end

%% Display holograms all
% close all
ax = [];
for tt = 1:size(H,3)
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    tiledlayout('flow')
    nexttile; imagesc(H(:,:,tt),[0,255]); axis image; colormap gray
    title(['Hologram ',num2str(tt)]); axis off; ax = [ax, gca];
    pause(0.1)
end
linkaxes(ax)

%% Reconstruct all holograms
for tt = 1:size(H,3)
    Y(:,:,tt) = AS_propagate_p(H(:,:,tt),Z(tt),n0,lambda,dx);
end

%% Display reconstructed phases (2 ROIs)
% close all
for tt = 1:size(H,3)
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    tiledlayout('flow')
    nexttile; imagesc(angle(Y(y1,x1,tt)),[-pi,pi]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 1']); axis off
    nexttile; imagesc(angle(Y(y2,x2,tt)),[-pi,pi]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 2']); axis off
    pause(0.1)
end

%% Invert phase (phase features should be positive) and normalize in 0-1 range
% For different cases may work better without inverting phases or for
% amplitude images or even for binarized amplitude/phase
for tt = 1:size(H,3)
    I(:,:,tt) = (-angle(Y(:,:,tt))+pi)/2/pi;
end

%% Find and apply affine transform to all holograms to correct their xy shift and scale
[I2,H2,tforms] = AutoAffineTransform(I,H);

%% Display corrected phases (2 ROIs)
% close all
for tt = 1:size(H,3)
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    tiledlayout('flow')
    nexttile; imagesc(-I2(y1,x1,tt),[-1,0]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 1']); axis off
    nexttile; imagesc(-I2(y2,x2,tt),[-1,0]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 2']); axis off
    pause(0.1)
end

%% Display corrected holograms (2 ROIs)
% close all
for tt = 1:size(H,3)
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    tiledlayout('flow')
    nexttile; imagesc(H2(y1,x1,tt),[0,255]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 1']); axis off
    nexttile; imagesc(H2(y2,x2,tt),[0,255]); axis image; colormap gray
    title(['Hologram ',num2str(tt),'; ROI 2']); axis off
    pause(0.1)
end

%% New propagation distances (because of rescalling the holograms)
Z2 = Z;
for tt = 1:size(H,3)-1
    tf = tforms{tt};
    Z2(tt) = Z(tt).*tf.T(1,1).^2;
end

%% Gercgberg Saxton phase retrieval
N_iter = 5; % Number of iterations (5 recommended)
useGPU = 1; % 1 if want to speed up calculations wth GPU (matlab need to suport used GPU)
Ygs = GSalgorithm(H2,Z2,n0,lambda,dx,N_iter,useGPU);

% Angular spectrum propagation for comparison
Yas = AS_propagate_p(H2(:,:,1),Z2(1),n0,lambda,dx);

%% Displaying retrieved phases
% close all
ax = [];
figure('Units','normalized','OuterPosition',[0,0,1,1])
imagesc(angle(Yas),[-pi,pi]); axis image; colormap gray
title('Angular spectrum phase'); ax = [ax,gca];
figure('Units','normalized','OuterPosition',[0,0,1,1])
imagesc(angle(Ygs),[-pi,pi]); axis image; colormap gray
title('GS AT phase'); ax = [ax,gca];
linkaxes(ax)

%% Creating gif file for readme
figure('Units','normalized','OuterPosition',[0,0,0.594,1])
tiledlayout(2,2,'TileSpacing','none','Padding','tight')
rng = [0,150];
for tt = [1:6,6:-1:1]
    nexttile(1); imagesc(H(y1,x1,tt),[0,255]); axis image; colormap gray
    text(10,10,['Hologram ',num2str(tt),'; ROI 1'],'FontSize',20,'Color','y')
    axis off;
    nexttile(2); imagesc(H(y2,x2,tt),[0,255]); axis image; colormap gray
    text(10,10,['Hologram ',num2str(tt),'; ROI 2'],'FontSize',20,'Color','y')
    axis off;
    nexttile(3); imagesc(H2(y1,x1,tt),[0,255]); axis image; colormap gray
    text(10,10,['Hologram ',num2str(tt),' after AT; ROI 1'],'FontSize',20,'Color','y')
    axis off;
    nexttile(4); imagesc(H2(y2,x2,tt),[0,255]); axis image; colormap gray
    text(10,10,['Hologram ',num2str(tt),' after AT; ROI 2'],'FontSize',20,'Color','y')
    axis off;
    pause(0.1)
    exportgraphics(gcf,"vid.gif","Append",true)
end

%% Saving phase images (ROI 1)
img = (angle(Yas(y1,x1))+pi)/2/pi;
imwrite(img,'AS_phase.png');
img = (angle(Ygs(y1,x1))+pi)/2/pi;
imwrite(img,'GS_AT_phase.png');