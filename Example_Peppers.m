% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
%% Code for reproducing results in Figure 4
p = 200; % Change to p=2 for Laplacian results
% Load Data
close all
addpath('MexNonlocal')
image=imread('peppers_masked.png');
U_sq=(image(:,:,1)==255).*(image(:,:,2)==255).*(image(:,:,3)==255);
image=double(image)./255;
[m,n,d]=size(image);
U_sq=logical(U_sq);
demim=image;
demim(repmat(U_sq,1,1,3))=1;

%% Compute initial inpainting
im_r=image(:,:,1);
im_g=image(:,:,2);
im_b=image(:,:,3);
known_pixels=[im_r(~U_sq),im_g(~U_sq),im_b(~U_sq)];
mu=mean(known_pixels);
Sigma=cov(known_pixels)+1e-8*eye(d);
R=chol(Sigma);
noise_pixels=normrnd(0,1,[3,sum(sum(U_sq))]);
noise_pixels=R'*noise_pixels;
noise_pixels=repmat(mu',1,size(noise_pixels,2))+noise_pixels;
noise_r=zeros(size(image,1),size(image,2));
noise_g=zeros(size(image,1),size(image,2));
noise_b=zeros(size(image,1),size(image,2));
noise_r(U_sq)=noise_pixels(1,:);
noise_g(U_sq)=noise_pixels(2,:);
noise_b(U_sq)=noise_pixels(3,:);
noise=noise_r;
noise(:,:,2)=noise_g;
noise(:,:,3)=noise_b;
im_out=image;
im_out(repmat(U_sq,1,1,3))=noise(repmat(U_sq,1,1,3));
% Save and show initial inpainting
imwrite(im_out,'peppers_init.png');
imagesc(im_out);
title('Initial image')
axis image;
axis off;

%% Compute solution using tight extension
input=reshape(im_out(:,:,1),[],1);
input(:,2)=reshape(im_out(:,:,2),[],1);
input(:,3)=reshape(im_out(:,:,3),[],1);
out=iterative_nonlocal_random_mask(input,U_sq(:),m,n,'plot',0,'patch_p',11,...
    'patch_r',30,'epsilon',1e-8,'nNeighbors',40,'min_Steps',11, 'p_max',p);
im_out=reshape(out(:,1),m,n);
im_out(:,:,2)=reshape(out(:,2),m,n);
im_out(:,:,3)=reshape(out(:,3),m,n);
% Save and show result
imwrite(im_out,'peppers_result.png');
figure,
imagesc(im_out);
title(['Result using tight extension with p=',num2str(p)])
axis image;
axis off;