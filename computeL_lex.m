% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Computes the vector of Llex values.
% INPUT:
% graph - adjacency matrix of the graph (size NxN, where N is the number
%         of the vertices) 
% U     - logical if vertex i is in the boundary (size Nx1)
% erg   - vertex values produced by ADMM(size N*nx1, where n is the
%           dimension of each value at one vertex)
% OUTPUT:
% S_max   - Llex vector
function S_max=computeL_lex(graph,U,erg)
% Determine non boundary vertices
free = 1:length(U);
free = free(~U);
% Random order for visiting vertices
per = randperm(length(free));
S_max = zeros(1,length(free));
for m=1:length(free)
    % Break paramater for exiting after the optimal vertex value is
    % determined
    i = free(per(m));
    % Adjacent vertex values
    adj_ind = (graph(:,i)~= 0);
    adj_weights = graph(adj_ind,i);
    tmp = erg(adj_ind,:);
    len = size(tmp,1);
    % Old value of S
    new_val = sqrt(max(sum((tmp-repmat(erg(i,:),len,1)).^2,2)...
        .*(adj_weights.^2)));
    if ~isempty(new_val)
        S_max(m) = new_val;
    end
end
% Sort Llex vector
S_max = sort(S_max,'descend')';
end
