A = importdata('exactsolution.txt');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);

tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);
colorbar

P = peaks(40);
%C = del2(P);
surf(P,C)
title('surf(P,del2(P))')