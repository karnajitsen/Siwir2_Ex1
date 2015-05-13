A= importdata('solution.txt');
M = A(:,1);
N = A(:,2);
P = A(:,3);
tri2 = delaunay(M,N);
trimesh(tri2,M,N,P);