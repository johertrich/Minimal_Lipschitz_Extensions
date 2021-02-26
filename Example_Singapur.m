% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
%% Code for reproducing results in Figure 6
p = 200; % Change to p=2 for Laplacian results
% Load Data
close all
image=imread('Singapur_masked.png');
U_sq=(image(:,:,1)==255).*(image(:,:,2)==0).*(image(:,:,3)==0);
image=double(image)./255;
[m,n,d]=size(image);
U_sq=logical(U_sq);
demim=image;
demim(repmat(U_sq,1,1,3))=1;
imwrite(demim,'withmask.png');

%% Compute initilization
out=iterative_nonlocal(image,U_sq,'plot',1,'sigma',.045,'patch_p',...
    15,'patch_r',45,'epsilon',1e-7,'nNeighbors',45,'p_max',p,'mode',1);
im_out=reshape(out(:,1),m,n);
im_out(:,:,2)=reshape(out(:,2),m,n);
im_out(:,:,3)=reshape(out(:,3),m,n);
% Save and show initial result
imwrite(im_out,'Singapur_init.png');
imagesc(im_out);
title('Initial image')
axis image;
axis off;

%% Compute solution
out=iterative_nonlocal(im_out,U_sq,'plot',1,'sigma',.045,'patch_p',15,...
    'patch_r',45,'epsilon',1e-7,'nNeighbors',45,'p_max',p);
im_out=reshape(out(:,1),m,n);
im_out(:,:,2)=reshape(out(:,2),m,n);
im_out(:,:,3)=reshape(out(:,3),m,n);
% Save and show final result
imwrite(im_out,'Singapur_result.png');
figure,
imagesc(im_out);
title(['Result using tight extension with p=',num2str(p)])
axis image;
axis off;
