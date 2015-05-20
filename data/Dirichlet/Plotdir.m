A = importdata('exactsolution_h_32.txt');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);

tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);
xlabel('X');
ylabel('Y');
title 'Exact Solution for Dirichlet Problem';

A = importdata('solution_h_32.txt');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);
figure
tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);
xlabel('X');
ylabel('Y');
title 'Approximated Solution for Dirichlet Problem';

