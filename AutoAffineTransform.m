function [OI,OH,tforms] = AutoAffineTransform(I,H,dp)
% Function that xy shifts and rescales images I(:,:,1:end-1) to match the
% I(:,:,end) image
% Inputs:
%   I - 3D matrix containing reconstructed holograms - normalized in 0-1
%       range 
%   H - input holograms (will be shifted according to affine transforms 
%       found for I)
%   dp - if > 0 - show matching results
% Outputs:
%   OI - xy shifted and and rescaled reconstructions to match the 
%        I(:,:,end)
%   OH - xy shifted and and rescaled holograms to match the I(:,:,end)
%   tforms - found affine transforms
% 
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
%
% Last modified: 05.01.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optional variables
thr = 1000; % threshold for detectSURFFeatures - higher threshold will 
            % result in smaller number of points
R = 200; % Radius (in pixels) in which the same features in different 
         % images may be present 

% Original image (will not be modified)
original = I(:,:,end);
% Detect features in original
ptsOriginal = detectSURFFeatures(original,'MetricThreshold',thr);
% Extract feature descriptors.
[featuresOriginal,validPtsOriginal] = extractFeatures(original,ptsOriginal);

for tt = 1:size(I,3)-1
    % Distorted image - will be xy shofted and rescaled to match original
    % image
    distorted = I(:,:,tt);
    % Detect features in distorted
    ptsDistorted = detectSURFFeatures(distorted,'MetricThreshold',thr);
    % Extract feature descriptors
    [featuresDistorted,validPtsDistorted] = extractFeatures(distorted,ptsDistorted);
    % Match features by using their descriptors.
    % indexPairs = matchFeatures(featuresOriginal,featuresDistorted);
    indexPairs = matchFeaturesInRadius(featuresOriginal,featuresDistorted,ptsDistorted.Location,ptsOriginal.Location,R);
    % Retrieve locations of corresponding points for each image
    matchedOriginal = validPtsOriginal(indexPairs(:,1));
    matchedDistorted = validPtsDistorted(indexPairs(:,2));
    
    % % Show putative point matches.
    % figure;
    % showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
    % title('Putatively matched points (including outliers)');
    
    % 
    [tform, inlierIdx] = estimateGeometricTransform2D(matchedDistorted,matchedOriginal,'similarity');
    
    inlierDistorted = matchedDistorted(inlierIdx,:);
    inlierOriginal = matchedOriginal(inlierIdx,:);
    if dp > 0
        % Display matching point pairs used in the computation of the transformation
        figure;
        showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted);
        title(['Matching points between image no. ', num2str(tt),...
            ' and image no. ',  num2str(size(I,3))]);
        legend('Points in last image', 'Points in n-th image');
    end

    Roriginal = imref2d(size(original));
    OI(:,:,tt) = imwarp(distorted,tform,'OutputView',Roriginal);
    OH(:,:,tt) = imwarp(H(:,:,tt),tform,'OutputView',Roriginal);
    tforms{tt} = tform;
end
OI(:,:,end+1) = I(:,:,end);
OH(:,:,end+1) = H(:,:,end);
end

