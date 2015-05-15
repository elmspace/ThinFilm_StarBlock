% Program IsoSurface, run by typing run("./main.m")


% Clear all previous variables
clear all;
% Getting the coordiante variables from the data, and making them
% into the correct format, mesh
A=importdata('./xyz.dat');
%getting the concentration from the data file
B=importdata('./ABCD.dat');

[X,Y,Z]=meshgrid(A(:,2),A(:,1),A(:,3));


% The size is used for the for loop
x_size=size(A(:,2),1);
y_size=size(A(:,1),1);
z_size=size(A(:,3),1);

% Dividing up the concentration values
VA=zeros(x_size,y_size,z_size);
VB=zeros(x_size,y_size,z_size);
VHA=zeros(x_size,y_size,z_size);
VHS=zeros(x_size,y_size,z_size);
% Taking the concentration values and putting them into the correct
% format, mesh-format
ii=1;
for i=1:x_size,
    for j=1:y_size,
        for k=1:z_size,
            VA(i,j,k)=B(ii,1);
            VB(i,j,k)=B(ii,2);
            VHA(i,j,k)=B(ii,3);
            VHS(i,j,k)=B(ii,4);
            ii=ii+1;
        end
    end
end
% clearing useless variables
clear i;
clear j;
clear k;
clear ii;
% Thres3old for the isosurf, 
cutA=0.35;
cutB=0.65;
cutHA=0.12;
cutHS=0.15;


axis vis3d;
view([1,0.3,0.2]);
axis off;
daspect('mode');

% B
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pB = patch(isosurface(X,Y,Z,VB,cutB),'FaceColor','green','EdgeColor','none');
qB = patch(isocaps(X,Y,Z,VB,cutB),'FaceColor','green','EdgeColor', ...
          'none');
alpha(pB,0.5);
alpha(qB,0.5);

% A    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pA = patch(isosurface(X,Y,Z,VA,cutA),'FaceColor','red','EdgeColor','none');
qA = patch(isocaps(X,Y,Z,VA,cutA),'FaceColor','red','EdgeColor', ...
         'none');
alpha(pA,0.7);
alpha(qA,0.7);

% HA    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pHA = patch(isosurface(X,Y,Z,VHA,cutHA),'FaceColor','blue','EdgeColor','none');
qHA = patch(isocaps(X,Y,Z,VHA,cutHA),'FaceColor','blue','EdgeColor', ...
           'none');
alpha(pHA,0.7);
alpha(qHA,0.7);

% HS    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pHS = patch(isosurface(X,Y,Z,VHS,cutHS),'FaceColor','yellow','EdgeColor','none');
qHS = patch(isocaps(X,Y,Z,VHS,cutHS),'FaceColor','yellow','EdgeColor', ...
           'none');
alpha(pHS,0.7);
alpha(qHS,0.7);
