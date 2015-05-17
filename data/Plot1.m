A = importdata('exactsolution_h_128.000000');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);

tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);

A = importdata('solution_h_128.000000');
X = A(:,1);
Y =A(:,2);
Z = A(:,3);
figure
tri = delaunay(X,Y);
trisurf(tri,X,Y,Z);
