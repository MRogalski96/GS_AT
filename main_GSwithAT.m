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

%% Input parameters
clear; close all; clc

% System parameters
Z = [11.4, 11.87, 12.4, 12.93, 13.45, 14]*1e3; % propagation distances [um]
lambda = 0.405; % wavelength [um]
dx = 2.4; % effective pixel size [um]
n0 = 1; % refractive index of medium (air)
nH = 6; % number of holograms

% GS parameters
N_iter = 5; % Number of iterations (5 recommended)
useGPU = 1; % 1 if want to speed up calculations wth GPU (matlab need to suport used GPU)

% Parameters for displaying results
dp = 1; % 0 - do not display results, 1 - display only important stuff, 2 - display all
x1 = 1001:1600; y1 = 851:1250; % ROI 1
x2 = 4001:4600; y2 = 1201:1600; % ROI 2

%% Code
% Load holograms
for tt = 1:nH
    H(:,:,tt) = double(imread(['./Data/H',num2str(tt),'.png']));
end


% Display loaded holograms (2 ROIs)
if dp > 0
    for tt = 1:nH
        figure('Units','normalized','OuterPosition',[0,0,1,1])
        tiledlayout('flow')
        nexttile; imagesc(H(y1,x1,tt),[0,255]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 1']); axis off
        nexttile; imagesc(H(y2,x2,tt),[0,255]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 2']); axis off
        pause(0.1)
    end
end

% Display holograms (full size)
if dp > 1
    ax = [];
    for tt = 1:nH
        figure('Units','normalized','OuterPosition',[0,0,1,1])
        tiledlayout('flow')
        nexttile; imagesc(H(:,:,tt),[0,255]); axis image; colormap gray
        title(['Hologram ',num2str(tt)]); axis off; ax = [ax, gca];
        pause(0.1)
    end
    linkaxes(ax)
end

% Reconstruct all holograms with AS method
for tt = 1:size(H,3)
    Y(:,:,tt) = AS_propagate_p(H(:,:,tt),Z(tt),n0,lambda,dx);
end

% Display reconstructed phases (2 ROIs)
if dp > 1
    for tt = 1:nH
        figure('Units','normalized','OuterPosition',[0,0,1,1])
        tiledlayout('flow')
        nexttile; imagesc(angle(Y(y1,x1,tt)),[-pi,pi]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 1']); axis off
        nexttile; imagesc(angle(Y(y2,x2,tt)),[-pi,pi]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 2']); axis off
        pause(0.1)
    end
end

% Invert phase (phase features should be positive) and normalize in 0-1 range
% For different cases may work better without inverting phases or for
% amplitude images or even for binarized amplitude/phase
for tt = 1:nH
    I(:,:,tt) = (-angle(Y(:,:,tt))+pi)/2/pi;
end

% Find and apply affine transform to all holograms to correct their xy shift and scale
[I2,H2,tforms] = AutoAffineTransform(I,H,dp-1);

% Display corrected phases (2 ROIs)
if dp > 1
    for tt = 1:nH
        figure('Units','normalized','OuterPosition',[0,0,1,1])
        tiledlayout('flow')
        nexttile; imagesc(-I2(y1,x1,tt),[-1,0]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 1']); axis off
        nexttile; imagesc(-I2(y2,x2,tt),[-1,0]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 2']); axis off
        pause(0.1)
    end
end

% Display corrected holograms (2 ROIs)
if dp > 0
    for tt = 1:nH
        figure('Units','normalized','OuterPosition',[0,0,1,1])
        tiledlayout('flow')
        nexttile; imagesc(H2(y1,x1,tt),[0,255]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 1']); axis off
        nexttile; imagesc(H2(y2,x2,tt),[0,255]); axis image; colormap gray
        title(['Hologram ',num2str(tt),'; ROI 2']); axis off
        pause(0.1)
    end
end

% New propagation distances (because of rescalling the holograms)
Z2 = Z;
for tt = 1:nH-1
    tf = tforms{tt};
    Z2(tt) = Z(tt).*tf.T(1,1).^2;
end

% Gercgberg Saxton phase retrieval
Ygs = GSalgorithm(H2,Z2,n0,lambda,dx,N_iter,useGPU);
% Angular spectrum propagation for comparison
Yas = AS_propagate_p(H2(:,:,1),Z2(1),n0,lambda,dx);

% Displaying retrieved phases
if dp > 0
    ax = [];
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    imagesc(angle(Yas),[-pi,pi]); axis image; colormap gray
    title('Angular spectrum phase'); ax = [ax,gca];
    figure('Units','normalized','OuterPosition',[0,0,1,1])
    imagesc(angle(Ygs),[-pi,pi]); axis image; colormap gray
    title('GS AT phase'); ax = [ax,gca];
    linkaxes(ax)
end