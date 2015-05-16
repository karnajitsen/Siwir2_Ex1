A = importdata('solution.txt');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);

tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);

