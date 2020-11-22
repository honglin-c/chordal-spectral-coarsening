function [U1,U2,D] = compute_diffusion_tensor(vertices,faces,opts)
% COMPUTE_DIFFUSION_TENSOR computes an averaged diffusion tensor on each
% face according to a certain \alpha parameter
% INPUT :
% vertices is a (nb_vertices)*3 matrix containing the coordinates of the
%   points on the shape
% faces is a (nb_faces)*3 matrix containing the indices of the vertices on
%   each face
% opts is a struct containing some fields:
%   opts.alpha is the anisotropic diffusion parameter
%   opts.curvature_smoothing is the parameter used to smooth the curvature
%   opts.use_sign asks whether to take into account the sign of curvature
%       or not (boolean). 
%   opts.normalize_tensor is a boolean to enforce the tensor to have a
%       determinant of 1
% OUTPUT:
% U1 is the (interpolated) minimal principal curvature direction on each
% triangle (vector of size (nb_faces)*3 )
% U2 is the (interpolated) maximal principal curvature direction on each
% triangle (vector of size (nb_faces)*3 )
% D is the (compressed) tensor on each triangle: it has size (nb_faces)*2,
% and D(i,1) represents the component of D in U1(i,:), D(i,2) represents
% the component of D in U2(i,:). The diffusion tensor is supposed diagonal
% REMARK :
% Note that due to the Finite Element Method, the tensor which is returned
% is not the tensor D but R^{T}DR, where R = [0,-1; 1,0] i.e. diagonal 
% elements are inverted

%% Default parameters
default_use_sign = 1;
default_normalize_tensor = 0;
default_alpha = 0;
default_curvature_smoothing = 10;
%% Check the input
if (nargin < 3)
   opts = [];
   opts.use_sign = default_use_sign;
   opts.normalize_tensor = default_normalize_tensor;
   opts.curvature_smoothing = default_curvature_smoothing;
   opts.alpha = 0;
else
    if ~isfield(opts,'alpha')
        opts.alpha = default_alpha;
    end
    if ~isfield(opts,'use_sign')
        opts.use_sign = default_use_sign;
    end
    if ~isfield(opts,'normalize_tensor')
        opts.normalize_tensor = default_normalize_tensor;
    end
    if ~isfield(opts,'curvature_smoothing')
        opts.curvature_smoothing = default_curvature_smoothing;
    end
end

%% Function used
% You can play with it to obtain different results
    function y = psi(x,coef,use_sign)
        if ~(use_sign)
            newx = x.^2;
        else
            newx = x;
        end
        y = exp(-repmat(coef,size(newx)).*newx);
    end
%% Compute the curvature
sub_options.curvature_smoothing = opts.curvature_smoothing;
sub_options.verb=0;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,normals] = compute_curvature(...
            transpose(vertices),transpose(faces),sub_options);
        clear Cmean Cgauss;
Umin = Umin';
Umax = Umax';
normals = normals';
%% Compute the average eigenvectors
[V1,V2] = interpolate_basis(Umin,Umax,normals,faces);
%build an orthonormal base for each face
edge1= vertices(faces(:,2),:)-vertices(faces(:,1),:);
edge2= vertices(faces(:,3),:)-vertices(faces(:,1),:);
edge1 = edge1./repmat(sqrt(sum(edge1.^2,2)),[1 3]);
edge2 = edge2./repmat(sqrt(sum(edge2.^2,2)),[1 3]);
edge2 = edge2 - repmat(sum(edge1.*edge2,2),[1 3]).*edge1;
%project [V1,V2] upon this base and orthonormalize them
U1=repmat(sum(V1.*edge1,2),[1 3]).*edge1 + ...
    repmat(sum(V1.*edge2,2),[1 3]).*edge2;
U2=repmat(sum(V2.*edge1,2),[1 3]).*edge1 + ...
    repmat(sum(V2.*edge2,2),[1 3]).*edge2;
U1=U1./repmat(sqrt(sum(U1.^2,2)),[1 3]);
U2=U2-repmat(sum(U2.*U1,2),[1 3]).*U1;
U2=U2./repmat(sqrt(sum(U2.^2,2)),[1 3]);
%% Compute the average eigenvalues of D
Cminmean=(1/3)*(Cmin(faces(:,1))+Cmin(faces(:,2))+Cmin(faces(:,3)));
Cmaxmean=(1/3)*(Cmax(faces(:,1))+Cmax(faces(:,2))+Cmax(faces(:,3)));
D = zeros(size(faces,1),2);
% compute the values (invert due to FEM derivation)
D(:,2)=psi(abs(Cminmean),opts.alpha,opts.use_sign);
D(:,1)=psi(abs(Cmaxmean),opts.alpha,opts.use_sign);

% renormalize the determinant if you need to
if(opts.normalize_tensor)
    %compute the determinant
    myDet = D(:,1).*D(:,2);
    % divide each coefficient by the square root of the determinant, so 
    %that the final determinant is 1
    for i = 1:2
        D(:,i) = D(:,i)./sqrt(myDet);
    end
end
end