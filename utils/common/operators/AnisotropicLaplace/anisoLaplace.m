function L = anisoLaplace(vertices,faces,alpha, options)
%  Anisotropic Laplace-Beltrami Spectrum
%   INPUT:
%   vertices is (number of vertices) x 3 matrix
%   faces is a (number of faces) x 3 matrix
%   options specifies different options :
%       n_eigenvalues : number of eigenvalues computed
%       verb : display progression in the computing
%       alpha : extension of anisotropy used (alpha >=0)
%       curvature_smoothing: size of the ring used to smooth curvature
%       normalize_areas: normalize or not the Voronoi areas
%       normalize_tensor: enforce (or not) det(tensor)=1
%       use_sign: take into account the sign of curvature (or not)
%       
%   OUTPUT:
%   PHI is a vector containing the ALB eigenfunctions, of size
%       (number of vertices)*(number of eigenvalues)
%   E is the vector of LB eigenvalues (by default of size 100 x 1)
%   Am is the (number of vertices)^2 diagonal matrix containing the
%   area of Voronoi cells around each vertex
%% default parameters
default_num_eigenvalues = 100; % number of eigenvalues computed
default_curvature_smoothing = 10;%size of the ring used to smooth curvature
default_alpha = 0; % value of alpha (anisotropic factor)
default_normalize_areas = 1; % normalize the Voronoi areas (enforce sum=1)
default_normalize_tensor = 0; %normalize anisotropic tensor (enforce det=1)
default_use_sign = 1; %use or don't use the sign of the curvature
default_verb = 1; % gives a written output in the console or not
%% parameters
if(nargin==4)
    if ~isfield(options,'verb')
        options.verb=default_verb;
    end
    if ~isfield(options,'n_eigenvalues')
        options.n_eigenvalues=default_num_eigenvalues;
    end
    if ~isfield(options,'curvature_smoothing')
        options.curvature_smoothing = default_curvature_smoothing;
    end
    if ~isfield(options,'alpha')
        options.alpha = default_alpha;
    end
    if ~isfield(options,'normalize_areas')
        options.normalize_areas = default_normalize_areas;
    end
    if ~isfield(options,'normalize_tensor')
        options.normalize_tensor = default_normalize_tensor;
    end
    if ~isfield(options,'use_sign')
        options.use_sign = default_use_sign;
    end
else
        options.n_eigenvalues=default_num_eigenvalues;
        options.verb=default_verb;
        options.curvature_smoothing = default_curvature_smoothing;
%         options.alpha = default_alpha;
        options.normalize_areas = default_normalize_areas;
        options.normalize_tensor = default_normalize_tensor;
        options.use_sign = default_use_sign;
end


%% basic quantities

num_vertices = size(vertices,1);
num_faces = size(faces,1);
options.alpha = alpha;
% alpha = options.alpha;
n_eigenvalues = options.n_eigenvalues;

%% detect boundary vertices

% Calculate the (directed) adjacency matrix. adjacency_matrix(m,n) = 1 if the oriented
% boundary of a triangle contains the directed edge from vertex m to vertex
% n, and 0 otherwise. This matrix is not quite symmetric due to boundary edges.
adjacency_matrix = sparse([faces(:,1); faces(:,2); faces(:,3)], ...
                         [faces(:,2); faces(:,3); faces(:,1)], ...
    	                 ones(3 * num_faces, 1), ...
                         num_vertices, num_vertices, 3 * num_faces);
if any(any(adjacency_matrix > 1))
    options.method = 'slow';
    faces = transpose(perform_faces_reorientation(vertices,faces,options));
     % error('Triangles must be oriented consistently.')
end

clear adjacency_matrix

%% compute Curvature-aware Laplace-Beltrami matrix
% if options.verb
% fprintf('Computing Anisotropic Laplace-Beltrami operator...');
% end

%% first compute inner face angles

angles = 0*faces; 

for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    pp = vertices(faces(:,i2),:) - vertices(faces(:,i1),:);
    qq = vertices(faces(:,i3),:) - vertices(faces(:,i1),:);
    % normalize the vectors
    pp = pp ./ repmat( max(sqrt(sum(pp.^2,2)),eps), [1 3] );
    qq = qq ./ repmat( max(sqrt(sum(qq.^2,2)),eps), [1 3] );
    % compute angles
    angles(:,i1) = acos(sum(pp.*qq,2));
end

%% Compute the anisotropic operator on each face.
% U1 is the direction of minimal curvature, U2 maximal, and D represents
% the values of the anisotropic tensor in this basis (it is actually
% transformed 
[U1,U2,D]=compute_diffusion_tensor(vertices,faces,options);
%% Build the ALB operator 
L = sparse(num_vertices,num_vertices);
% non diagonal terms
for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    % we only focus on (i1,i1) and on (i1,i2)
    e1 = vertices(faces(:,i3),:)-vertices(faces(:,i2),:);
    e2 = vertices(faces(:,i1),:)-vertices(faces(:,i3),:);
    % normalize e1 and e2
    e1 = e1./repmat(sqrt(sum(e1.^2,2)),[1 3]);
    e2 = e2./repmat(sqrt(sum(e2.^2,2)),[1 3]);
    %edge factor for matrix L (non-diagonal factor (i1,i2) ):
    % corresponds to -(1/2)<De1,e2>/sin(alpha_{12})
    factore=-(1/2)*(D(:,1).*(sum(e1.*U1(:,:),2)).*(sum(e2.*U1(:,:),2)) +...
                D(:,2).*(sum(e1.*U2(:,:),2)).*(sum(e2.*U2(:,:),2)))...
                ./sin(angles(:,i3));
    %diagonal factor (i1,i1):
    factord=-(1/2)*(D(:,1).*(sum(e1.*U1(:,:),2).^2) + ...
        D(:,2).*(sum(e1.*U2(:,:),2).^2)).*(cot(angles(:,i2))+cot(angles(:,i3)));
    L = L + sparse([faces(:,i1); faces(:,i2); faces(:,i1)],...
        [faces(:,i2) ;faces(:,i1) ;faces(:,i1)],[factore ;factore ;factord],...
        num_vertices, num_vertices);
end
clear angles;