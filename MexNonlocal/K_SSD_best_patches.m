function out = K_SSD_best_patches(input, varargin)
% [Samples] = K_SSD_best_patches(Img, ...)
% 'SampleSize'
% 'PatchWidth'
% 'SearchWindowRadius'
% 'Mask'
% 

size_of_the_sample= single(40);
patch_width=single(3);
search_window_radius=single(15);


p = inputParser;
addRequired(p, 'input', @isnumeric);
addOptional(p, 'SampleSize', size_of_the_sample, @isnumeric);
addOptional(p, 'PatchWidth', patch_width, @isnumeric);
addOptional(p, 'SearchWindowRadius', search_window_radius, @isnumeric);
addOptional(p, 'Mask', @isnumeric);

parse(p, input, varargin{:});

initialtype = class(input);
input=single(input);

size_of_the_sample= single(p.Results.SampleSize);
patch_width= single(p.Results.PatchWidth);
search_window_radius= single(p.Results.SearchWindowRadius);
if nargin()==9
    Mask= single(p.Results.Mask);
end

if ndims(input)==3 && size(input,3)==3
    if nargin==7
        out = SSD_similarity_3D(input, size_of_the_sample, patch_width,...
            search_window_radius);
    else
        out = SSD_similarity_3DM(input, size_of_the_sample, patch_width,...
            search_window_radius,Mask);
    end
else
    if nargin==7
        out = SSD_similarity(input, size_of_the_sample, patch_width,...
            search_window_radius);
    else
        out = SSD_similarityM(input, size_of_the_sample, patch_width,...
            search_window_radius,Mask);
    end
end

out = cast(out, initialtype);
end
