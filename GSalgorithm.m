function [uout,uin] = GSalgorithm(H,z,n0,lambda,dx,iter,useGPU)
% Function that performs Gerchberg-Saxton hologram reconstruction from
% several digital in-line holograms collected with different parameters 
% defocus/propagation distances.
% 
% Inputs:
%     H - 3D matrix containing a stack of holograms
%     z - vector containing propagation distances for all holograms -
%         distance between hologram plane and object plane
%         z(1) for H(:,:,1), z(2) for H(:,:,2) etc.
%     n0 - free space refractive index (for air n0 = 1)
%     lambda - wavelength
%     dx - camera pixel size in the object plane (CamPixSize/magnification)
%     iter - number of iterations (5 is recommended)
%     useGPU - = 1 if want to use GPU optimization (may not work for older 
%              GPUs or for too large H)
% Outputs:
%     uout - reconstructed complex optical field at the object plane
%     uin - reconstructed complex optical field at the H(:,:,1) plane
%  
% Important: all 'length' variables (z, lambda, dx), should be provided in 
% the same units (e.g., um)
% 
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
%
% Last modified: 05.01.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert holograms to GPU
if useGPU == 1
    H = gpuArray(H);
end
% Amplitude = square root of the holograms
amp = sqrt(H);


% Precompute kernels for faster AS propagation
[Ny,Nx,Nh] = size(H);
dfx = 1/Nx/dx; fx = -Nx/2*dfx : dfx : (Nx/2-1)*dfx;
dfy = 1/Ny/dx; fy = -Ny/2*dfy : dfy : (Ny/2-1)*dfy;
kernel = zeros(Ny,Nx,Nh);
if useGPU == 1
    kernel = gpuArray(kernel);
end
zz(1:Nh+1) = 0;
for hh = 1:Nh+1
    % Precompute propagation distances
    if hh == 1
        % From last plane to first plane
        zz(hh) = z(end) - z(1);
    elseif hh == Nh+1
        % From first plane to object plane
        zz(hh) = z(1);
    else
        % From hh-1 plane to hh plane
        zz(hh) = z(hh-1) - z(hh);
    end
    k = 2*pi/lambda;
    p = fftshift(k*abs(zz(hh))*sqrt(n0^2 - lambda^2*...
        (ones(Ny,1)*(fx.^2)+(fy'.^2)*ones(1,Nx))));
    p = p - p(1,1); % Subtract to have phase oscilating around 0
    kernel(:,:,hh) = exp(1i*p); % Precomputed kernel
end


% GS algorithm
uin = amp(:,:,1);
for tt = 1:iter
    for ss = 2:size(H,3)
        % Propagate optical field from ss-1 to ss hologram plane
        uin = AS_propagate(uin,zz(ss),kernel(:,:,ss));
        % Apply amplitude constraint (replace amplitude with measured one)
        uin = uin./abs(uin).*amp(:,:,ss);
    end
    % Propagate optical field from last to first hologram plane
    uin = AS_propagate(uin,zz(1),kernel(:,:,1));
    % Apply amplitude constraint (replace amplitude with measured one)
    uin = uin./abs(uin).*amp(:,:,1);
end
% Final propagation to the object plane
uout = AS_propagate(uin,zz(end),kernel(:,:,end));

if useGPU == 1
    uout = gather(uout);
end

end

function uout = AS_propagate(uin,z,kernel)
% Angular spectrum propagation
if  z<0 
    ftu = kernel.*(fft2(conj(uin)));
    uout = conj(ifft2(ftu));
else
    ftu = kernel.*fft2(uin);
    uout = ifft2(ftu);
end
end