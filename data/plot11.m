%for i=0:6
A= importdata('solution.txt');
M = A(:,1);
N = A(:,2);
P = A(:,3);
tri2 = delaunay(M,N);
figure
trimesh(tri2,M,N,P);
title "Down";
%end
% B= importdata('exactsolution.txt');
% M = B(:,1);
% N = B(:,2);
% P = B(:,3);
% tri2 = delaunay(M,N);
% figure
% trimesh(tri2,M,N,P);
% title 'Exact Solution';