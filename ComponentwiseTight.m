% This code belongs to the paper
%
% M. Bačák, J. Hertrich, S. Neumayer and G. Steidl.
% Minimal Lipschitz and ∞-Harmonic Extensions of Vector-Valued Functions on Finite Graphs.
% Information and Inference: A Journal of the IMA, vol 9, pp. 935–959, 2020.
% 
% Please cite the paper, if you use this code.
%
% Minimizes the norm of the Laplace-Operator
% INPUT:
% g         - adjacency matrix of graph (size NxN, N number of vertices)
% f         - values on given vertices
% maske     - area on which the extension should be computed
% OPTIONAL:
% delta_t       - step size in diffusion equation
% epsilon       - tolerance
% maxIterations - maximum number of iterations
% OUTPUT:
% out           - extended function values
function out=ComponentwiseTight(g,f,maske,varargin)
p=inputParser;
addOptional(p,'delta_t',0.4);
addOptional(p,'epsilon',1e-5);
addOptional(p,'maxIterations',10000);

parse(p,varargin{:});
s=p.Results;
LapFun = @(g,f)infLaplacian(g,f,0.5);

% Restrict graph to relevant part
relevant=((g*maske+maske)~=0);
graph_relevant=g(relevant,relevant);
f_n=f(relevant,:);
maske_relevant=maske(relevant);
f_min=f_n;

if isa(s.delta_t,'function_handle')==1
    delta_t=s.delta_t(graph_relevant);
else
    delta_t=s.delta_t;
end

Lap=LapFun(graph_relevant,f_n);
infLap_norm=sum(abs(Lap),2);
iterations=0;
abbruchkriterium=max(abs(infLap_norm.*maske_relevant));
abbruch_min=abbruchkriterium;

% Iterate until tolerance is satisfied or maximum number of iterations
% is reached
while abbruchkriterium>s.epsilon && iterations<s.maxIterations
    % Euler updates
    f_n=f_n+delta_t*(repmat(maske_relevant,1,size(f_n,2)).*Lap);
    
    % Recompute Laplacian
    Lap=LapFun(graph_relevant,f_n);
    infLap_norm=sum(abs(Lap),2);
    iterations=iterations+1;
    abbruchkriterium=sum(abs(infLap_norm.*maske_relevant))/size(maske_relevant,1);
    if abbruchkriterium<abbruch_min
        abbruch_min=abbruchkriterium;
        f_min=f_n;
    end
    if mod(iterations,10000)==0
        disp(abbruchkriterium);
        disp(iterations);
    end
end
if iterations==s.maxIterations
    abbruchkriterium=abbruch_min;
    f_n=f_min;
end
disp(abbruchkriterium)
disp(iterations)
out=f;
out(relevant,:)=f_n;
end

function lap=infLaplacian(graph,f,alpha)
% Returns the infinity Laplacian L_{omega,inf} in each component as defined
% in Elmoataz, Toutain, Tenbrinck 2015
% INPUT:
% g     - adjacency matrix of a weighted graph
% f     - input function
% alpha - weight
% OUTPUT:
% lap   - L_{omega,inf} f
beta=1-alpha;
[n,d]=size(f);
lap=zeros(n,d);
wurzel_Gewichte=graph.^(0.5);
for i=1:d
    fi=f(:,i);
    graph_mal_f_u=diag(sparse(fi))*wurzel_Gewichte;
    graph_mal_f_v=wurzel_Gewichte*diag(sparse(fi));
    upwind_plus=max(max(graph_mal_f_v-graph_mal_f_u,0),[],2);
    upwind_minus=max(max(graph_mal_f_u-graph_mal_f_v,0),[],2);
    lap(:,i)=alpha.*upwind_plus-beta.*upwind_minus;
end
end
