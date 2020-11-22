function [V1,V2] = interpolate_basis(U1,U2,normals,faces)
% interpolate_basis  interpolates an orthonormal basis on each face, on a
% point given by its length ratio lratio (from the first vertex of the
% face) and its angular ratio aratio (from v1v2 to v1v3)

% INPUT : 
% U1 is a (nb_vertices)*3 matrix representing the 1st basis vector on each
% vertex (usually Umin for principal curvature directions)
% U2 is a (nb_vertices)*3 matrix representing the 2nd basis vector on each
% vertex (usually Umax for principal curvature directions)
% normals 
% faces is a (nb_faces)*3 matrix representing the faces of the shape

% OUTPUT
% V1 is a (nb_faces)*3 vector (the first vector of the basis) (usually 
% maximal curvature direction)
% V2 is a (nb_faces)*3 vector (the second vector of the basis) (usually 
% minimal curvature direction)


%% Make such that (U1,U2) is direct (U1^U2 is oriented along) and compute 3rd vector
U3 = cross(U1,U2,2);
d = sum(U3.*normals,2);
I = find(d<0);

if(numel(I)>0)
    U2(I,:) = -U2(I,:);
    U3 = cross(U1,U2,2);
end

%% Average the vectors
V1 = (1/3)*(U1(faces(:,1),:)+U1(faces(:,2),:)+U1(faces(:,3),:));
V2 = (1/3)*(U2(faces(:,1),:)+U2(faces(:,2),:)+U2(faces(:,3),:));
end

