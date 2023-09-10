clear; clc; close;
%%
load("dalekosaur\object.mat");
patch('vertices', Xo', 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'k');
hold on
axis equal;
xlabel('Xo-axis'); ylabel('Yo-axis'); zlabel('Zo-axis');
%%
im = imread("images\IMG_8299.jpg");
ObjectDirectory = 'dalekosaur';
[impoints, objpoints3D] = clickPoints(im, ObjectDirectory);

%% plotting correspondance points
figureIndex = 1;
figure(4)
figureIndex = figureIndex + 1;
imshow(im); 
hold on;
plot( impoints(:,1), impoints(:,2), 'b.');
hold off

figure(5)
figureIndex = figureIndex + 1;
patch('vertices', Xo', 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'k');
hold on;
axis vis3d;
axis equal;
plot3( objpoints3D(:,1), objpoints3D(:,2), objpoints3D(:,3), 'b.');
xlabel('Xo-axis'); ylabel('Yo-axis'); zlabel('Zo-axis');
hold off;

%% calculating M

M = estimateCameraProjectionMatrix(impoints,objpoints3D)

%% KRt
A = M(:,1:3);
b = M(:,4);
C = A*A';

lambda1 = 1/sqrt(C(3,3));
lambda2 = -1/sqrt(C(3,3));

lambdasq = 1/C(3,3);

Xc = C(1,3)/C(3,3);
Yc = C(2,3)/C(3,3);
fy = sqrt(lambdasq*C(2,2)-(Yc*Yc));
alpha = (1/fy)*(lambdasq*C(1,2) - (Xc*Yc));
fx = sqrt(lambdasq*C(1,1)-(alpha*alpha)-(Xc*Xc));

K = [fx, alpha, Xc; 0, fy, Yc; 0, 0, 1]

R1 = lambda1*inv(K)*A;
R2 = lambda2*inv(K)*A;
R = zeros(3,3);
lambda = 0;
if det(R1) > 0
    R = R1
    lambda = lambda1;
else
    R = R2
    lambda = lambda2;
end

t = lambda*inv(K)*b

%% Verifying Result

[points,w] = size(objpoints3D);
estPoints = [];
sumsq = 0;
for i = 1:points
    p = objpoints3D(i,:);
    p = [p';1];
    out = K*[R,t]*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe;
    ye = ype/wpe;
    
    estPoints = [estPoints; [xe,ye]];
    sumsq = ((impoints(i,1)-xe)*(impoints(i,1)-xe)) + ((impoints(i,2)-ye)*(impoints(i,2)-ye));
end

meansq = sumsq/points

figure(3)
figureIndex = figureIndex + 1;
imshow(im); 
hold on;
plot( impoints(:,1), impoints(:,2), 'b.');
plot(estPoints(:,1), estPoints(:,2), 'ro');
hold off

%% superimposed mesh

[w,points] = size(Xo);
xtheta = 0;
ytheta = 0;
ztheta = 0;
tx = 0;
ty = 0;
tz = 0;
rotx = [1,0,0,0;0,cosd(xtheta),-sind(xtheta),0;0,sind(xtheta),cosd(xtheta),0;0,0,0,1];
roty = [cosd(ytheta),0,sind(ytheta),0;0,1,0,0;-sind(ytheta),0,cosd(ytheta),0;0,0,0,1];
rotz = [cosd(ztheta),-sind(ztheta),0,0;sind(ztheta),cosd(ztheta),0,0;0,0,1,0;0,0,0,1];
trans = [1,0,0,tx;0,1,0,ty;0,0,1,tz;0,0,0,1];
newPoints = [];

for i = 1:points
    p = [Xo(:,i);1];
    p = roty*rotx*rotz*trans*p;
    %p = [p;1];

    out = K*[R,t]*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe;
    ye = ype/wpe;
    newPoints = [newPoints; [xe,ye]];
end



figure(figureIndex)
figureIndex = figureIndex + 1;
imshow(im);
hold on;
patch('vertices', newPoints, 'faces', Faces, 'facecolor', 'n', 'edgecolor', 'b');
hold off;
%%2
cameraCalibrator('checkerboard\',37)

%%

Kchecker = zeros(3,3);
fc = cameraParams.FocalLength;
cc = cameraParams.PrincipalPoint;
alpha_c = cameraParams.Skew;
Kchecker(1,1) = fc(1);
Kchecker(1,2) = alpha_c;
Kchecker(1,3) = cc(1);
Kchecker(2,2) = fc(2);
Kchecker(2,3) = cc(2);
Kchecker(3,3) = 1;

%% 3
figureIndex = 8;
galleryImage = imread("images\IMG_8210.jpg");
xtheta = 20;
ytheta = 80;
ztheta = 30;
tx = 10;
ty = 1000;
rotx = [1,0,0,0;0,cosd(xtheta),-sind(xtheta),0;0,sind(xtheta),cosd(xtheta),0;0,0,0,1];
roty = [cosd(ytheta),0,sind(ytheta),0;0,1,0,0;-sind(ytheta),0,cosd(ytheta),0;0,0,0,1];
rotz = [cosd(ztheta),-sind(ztheta),0,0;sind(ztheta),cosd(ztheta),0,0;0,0,1,0;0,0,0,1];
scale = [.000001,0,0;0,.000001,0;0,0,.000001];
rot = rotz*rotx*roty;
trans = [1,0,0,tx;0,1,0,ty;0,0,1,tz;0,0,0,1];
x_projected = [];
x_transformed = [];

for i = 1:points
    p = [Xo(:,i);1];
    p = [R,t]*rot*p;
    %p = [p;1];
    x_transformed = [x_transformed,p(1:3)];

end


for i = 1:points
    p = [x_transformed(:,i)];
    out = K*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(8)

imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'k');
displayLit(x_projected,x_transformed,Faces,[1,1,0],pointsInFront);
hold off

% Image 1 1 kchecker
x_projected = [];
tx = 200;
ty = 800;
for i = 1:points
    p = [x_transformed(:,i)];
    out = Kchecker*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(9)
figureIndex = figureIndex + 1;
imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'k');
displayLit(x_projected,x_transformed,Faces,[1,1,0],pointsInFront);
hold off

%% Image 1 2

galleryImage = imread("IMG_8210.jpg");
xtheta = 20;
ytheta = -80;
ztheta = 30;
tx = -1000;
ty = 1000;
rotx = [1,0,0,0;0,cosd(xtheta),-sind(xtheta),0;0,sind(xtheta),cosd(xtheta),0;0,0,0,1];
roty = [cosd(ytheta),0,sind(ytheta),0;0,1,0,0;-sind(ytheta),0,cosd(ytheta),0;0,0,0,1];
rotz = [cosd(ztheta),-sind(ztheta),0,0;sind(ztheta),cosd(ztheta),0,0;0,0,1,0;0,0,0,1];
scale = [.000001,0,0;0,.000001,0;0,0,.000001];
rot = rotz*rotx*roty;
trans = [1,0,0,tx;0,1,0,ty;0,0,1,tz;0,0,0,1];
x_projected = [];
x_transformed = [];

for i = 1:points
    p = [Xo(:,i);1];
    p = [R,t]*rot*p;
    %p = [p;1];
    x_transformed = [x_transformed,p(1:3)];

end


for i = 1:points
    p = [x_transformed(:,i)];
    out = K*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(10)

imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,0],pointsInFront);
hold off

%  Image 1 2 Kchecker
x_projected = [];
tx = -800;
ty = 800;
for i = 1:points
    p = [x_transformed(:,i)];
    out = Kchecker*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(11)

imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,0],pointsInFront);
hold off

%% Image 1 3

galleryImage = imread("IMG_8210.jpg");
xtheta = 20;
ytheta = 120;
ztheta = 30;
tx = 1000;
ty = 1000;
rotx = [1,0,0,0;0,cosd(xtheta),-sind(xtheta),0;0,sind(xtheta),cosd(xtheta),0;0,0,0,1];
roty = [cosd(ytheta),0,sind(ytheta),0;0,1,0,0;-sind(ytheta),0,cosd(ytheta),0;0,0,0,1];
rotz = [cosd(ztheta),-sind(ztheta),0,0;sind(ztheta),cosd(ztheta),0,0;0,0,1,0;0,0,0,1];
scale = [.000001,0,0;0,.000001,0;0,0,.000001];
rot = rotz*rotx*roty;
trans = [1,0,0,tx;0,1,0,ty;0,0,1,tz;0,0,0,1];
x_projected = [];
x_transformed = [];

for i = 1:points
    p = [Xo(:,i);1];
    p = [R,t]*rot*p;
    %p = [p;1];
    x_transformed = [x_transformed,p(1:3)];

end


for i = 1:points
    p = [x_transformed(:,i)];
    out = K*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(12)

imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[1,1,0],pointsInFront);
hold off

% Image 1 3 Kchecker
x_projected = [];
tx = 1500;
ty = 700;
scale = [0.001,0,0;0,0.001,0;0,0,.001];
for i = 1:points
    p = [x_transformed(:,i)];
    out = Kchecker*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(13)

imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[1,1,0],pointsInFront);
hold off

%% Image 2 1
galleryImage = imread("uganda.jpg");
xtheta = 20;
ytheta = 0;
ztheta = 30;
tx = 350;
ty = 175;
rotx = [1,0,0,0;0,cosd(xtheta),-sind(xtheta),0;0,sind(xtheta),cosd(xtheta),0;0,0,0,1];
roty = [cosd(ytheta),0,sind(ytheta),0;0,1,0,0;-sind(ytheta),0,cosd(ytheta),0;0,0,0,1];
rotz = [cosd(ztheta),-sind(ztheta),0,0;sind(ztheta),cosd(ztheta),0,0;0,0,1,0;0,0,0,1];
scale = [.000001,0,0;0,.000001,0;0,0,.000001];
rot = rotz*rotx*roty;
trans = [1,0,0,tx;0,1,0,ty;0,0,1,tz;0,0,0,1];
x_projected = [];
x_transformed = [];

for i = 1:points
    p = [Xo(:,i);1];
    p = [R,t]*rot*p;
    %p = [p;1];
    x_transformed = [x_transformed,p(1:3)];

end


for i = 1:points
    p = [x_transformed(:,i)];
    out = K*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(14)
imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,1],pointsInFront);
hold off

% image 2 1 kchecker 
x_projected = [];
tx = 750;
ty = 0;

for i = 1:points
    p = [x_transformed(:,i)];
    out = Kchecker*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(15)
imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,1],pointsInFront);
hold off

%% image 2 2 k
galleryImage = imread("uganda.jpg");
xtheta = 20;
ytheta = 0;
ztheta =0;
tx = 1400;
ty = 200;
rotx = [1,0,0,0;0,cosd(xtheta),-sind(xtheta),0;0,sind(xtheta),cosd(xtheta),0;0,0,0,1];
roty = [cosd(ytheta),0,sind(ytheta),0;0,1,0,0;-sind(ytheta),0,cosd(ytheta),0;0,0,0,1];
rotz = [cosd(ztheta),-sind(ztheta),0,0;sind(ztheta),cosd(ztheta),0,0;0,0,1,0;0,0,0,1];
rot = rotz*rotx*roty;
trans = [1,0,0,tx;0,1,0,ty;0,0,1,tz;0,0,0,1];
x_projected = [];
x_transformed = [];

for i = 1:points
    p = [Xo(:,i);1];
    p = [R,t]*rot*p;
    %p = [p;1];
    x_transformed = [x_transformed,p(1:3)];

end


for i = 1:points
    p = [x_transformed(:,i)];
    out = K*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(16)
imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,1],pointsInFront);
hold off

% image 22 k checker
x_projected = [];
tx = 1600;
ty = 0;

for i = 1:points
    p = [x_transformed(:,i)];
    out = Kchecker*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(17)
imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,1],pointsInFront);
hold off

%% image 2 3 k 
galleryImage = imread("uganda.jpg");
xtheta = 10;
ytheta = -90;
ztheta =20;
tx = 1400;
ty = 900;
rotx = [1,0,0,0;0,cosd(xtheta),-sind(xtheta),0;0,sind(xtheta),cosd(xtheta),0;0,0,0,1];
roty = [cosd(ytheta),0,sind(ytheta),0;0,1,0,0;-sind(ytheta),0,cosd(ytheta),0;0,0,0,1];
rotz = [cosd(ztheta),-sind(ztheta),0,0;sind(ztheta),cosd(ztheta),0,0;0,0,1,0;0,0,0,1];
scale = [.000001,0,0;0,.000001,0;0,0,.000001];
rot = rotz*rotx*roty;
trans = [1,0,0,tx;0,1,0,ty;0,0,1,tz;0,0,0,1];
x_projected = [];
x_transformed = [];

for i = 1:points
    p = [Xo(:,i);1];
    p = [R,t]*rot*p;
    %p = [p;1];
    x_transformed = [x_transformed,p(1:3)];

end


for i = 1:points
    p = [x_transformed(:,i)];
    out = K*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end

pointsInFront = isinfront(x_transformed,Faces);

figure(18)
imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,1],pointsInFront);
hold off

% image 2 3 kchecker
x_projected = [];
tx = 1800;
ty = 600;

for i = 1:points
    p = [x_transformed(:,i)];
    out = Kchecker*p;
    xpe = out(1,1);
    ype = out(2,1);
    wpe = out(3,1);

    xe = xpe/wpe + tx;
    ye = ype/wpe + ty;
    x_projected = [x_projected; [xe,ye]];
end
pointsInFront = isinfront(x_transformed,Faces);

figure(19)
imshow(galleryImage)
hold on
%patch('vertices', x_projected, 'faces', Faces, 'facecolor', 'w', 'edgecolor', 'b');
displayLit(x_projected,x_transformed,Faces,[0,1,1],pointsInFront);
hold off
