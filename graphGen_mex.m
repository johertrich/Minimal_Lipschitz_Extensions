% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Computes an adjacency matrix of a graph on an image based on patch similarity
% INPUT:
% img               - image
% mask              - unknown pixels
% OPTIONAL:
% patchSize         - A patch around a pixel (i,j) is defined by 
%                     [i-patchSize:i+patchSize,j-patchSize:j+patchSize]
% nNeighbors        - number of neighbors in the graph
% radius            - Search for similar patches in a neighborhood of size 2*radius+1
% OUTPUT:
% out               - adjacency matrix of the resulting graph
function out=graphGen_mex(img,mask,m,n,varargin)
    ip=inputParser;
    addOptional(ip,'patchSize',7);
    addOptional(ip,'nNeighbors',8);
    addOptional(ip,'radius',15);
    parse(ip,varargin{:});
    
    s=ip.Results;
    p=s.patchSize;
    r=s.radius;
    k=s.nNeighbors;
    
    [~,d]=size(img);
    img2 = zeros(m,n,d);
    for i =1:d
        img2(:,:,i) = reshape(img(:,i),m,n);
    end
    % Mask2 for finding Pixels with full Patch
    mask = reshape(mask,m,n);
    [indx,indy]=find(mask==1);
    nPixel=size(indx,1);
    if rem(p,2) ==0
        error('Only odd numbers allowed!')
    end   
    
    % Compute nearest neighbors with MexImplementation
    if d ==1 || d == 3
        img2 = img2*255;
    else
        error('Only implemented d=1,3')
    end
    neigh = K_SSD_best_patches(img2,'SampleSize',k,...
        'PatchWidth', p,'SearchWindowRadius', r);
    largestWeights = double(squeeze(neigh(:,:,:,1)));
    largestWeights = exp(-largestWeights./(largestWeights(:,:,20)+1e-9));
    %Create indices
    ind1 = squeeze(neigh(:,:,:,2))+1;
    ind2 = squeeze(neigh(:,:,:,3));
    ind1 = ind1+m*ind2;
    data = zeros(k*nPixel,3);
    counter = 0;
    % Construct graph
    for i=1:nPixel % Iterate over all unknown points
        x=indx(i);
        y=indy(i);
        % Create data for graph
        data(counter+1:counter+k,1)=squeeze(ind1(x,y,:));
        data(counter+1:counter+k,2)=(x + m*(y-1))*ones(k,1);
        data(counter+1:counter+k,3)=squeeze(largestWeights(x,y,:));
        counter = counter + k;
    end
    out=sparse(data(:,1),data(:,2),data(:,3),m*n,m*n);
    out= out + out';
end
