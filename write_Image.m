% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Save and plot a vectorized image
% INPUT:
%   filename - path for saving
%   values - vectorized image
%   m,n - size of the image
%   plot - 1 for plotting the image, 0 or no input else
function write_Image(filename,values,m,n,plot)
    rgb_vector=reshape(values,3,size(values,1)/3)';
    image=reshape(rgb_vector(:,1),m,n);
    image(:,:,2)=reshape(rgb_vector(:,2),m,n);
    image(:,:,3)=reshape(rgb_vector(:,3),m,n);
    imwrite(image,filename)
    if nargin==5&&plot==1
        figure;
        imagesc(image);
        axis equal
        axis off
    end
end
