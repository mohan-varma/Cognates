%% Code to find Cognates:
clc
close all
% Given a four bar mechanism O1-A-B-O2 with coupler point P. This code
% finds the cognates (two other four bar mechanisms and the corrsponding
% coupler points that will trace the exact same coupler curve.

% Inputs: Link lengths lg (ground), O1-A, A-B, O2-B, A-P and B-P and angle
% of link O1-A angA
% Outputs: lengths lg1, O11-A1, A1-B1, O21-B1, A1-P1 and B1-P1 and 
%                  lg2, O12-A2, A2-B2, O22-B2, A2-P2 and B2-P2 and
% Inputs:
lg = 10; O1_A = 4; A_B = 7; O2_B = 8; A_P = 4; B_P = 4; angA = 100;

% Get coordinates all points O1, A, B, P, O2 such that AB is horizontal.
Pos_O1 = [0;0;0];
Pos_A = O1_A*[cosd(angA);sind(angA);0];
Pos_O2 = [lg;0;0];

[Ix,Iy] = circcirc(Pos_A(1,1),Pos_A(2,1),A_B,Pos_O2(1,1),Pos_O2(2,1),O2_B);
if(Iy(1)>Iy(2))
    Pos_B  = [Ix(1);Iy(1);0];
else
    Pos_B  = [Ix(2);Iy(2);0];
end
[Ix,Iy] = circcirc(Pos_A(1,1),Pos_A(2,1),A_P,Pos_B(1,1),Pos_B(2,1),B_P);
if(Iy(1)>Iy(2))
    Pos_P  = [Ix(1);Iy(1);0];
else
    Pos_P  = [Ix(2);Iy(2);0];
end

pointsx = [Pos_O1(1,1);Pos_A(1,1);Pos_B(1,1);Pos_P(1,1);Pos_A(1,1);Pos_B(1,1);Pos_O2(1,1)];
pointsy = [Pos_O1(1,1);Pos_A(2,1);Pos_B(2,1);Pos_P(2,1);Pos_A(2,1);Pos_B(2,1);Pos_O2(2,1)];
figno = 0;
%Figure 1
figno = figno + 1;
figure(figno)
plot(pointsx,pointsy,'o')
hold on
line(pointsx,pointsy)
hold off
grid on
axis equal
title('Given Mechanism')


%% Step 1: move pivot points such that O1-A-B-O2 are collinear. Moved pivot 
% points are named with additional '_vert' indicating vertics of triangle 
% in the Cayley diagram. 
dir1 = (Pos_A-Pos_B)/norm(Pos_A-Pos_B);
Pos_O1_Cayley = Pos_A + O1_A*dir1;
dir2 = -dir1;
Pos_O2_Cayley = Pos_B + O2_B*dir2;
points3x = [Pos_O1_Cayley(1,1);Pos_A(1,1);Pos_B(1,1);Pos_O2_Cayley(1,1)];
points3y = [Pos_O1_Cayley(2,1);Pos_A(2,1);Pos_B(2,1);Pos_O2_Cayley(2,1)];
%Figure 2
figno = figno + 1;
figure(figno)
plot(pointsx,pointsy,'o')
hold on
plot(points3x,points3y,'ro')
line(pointsx,pointsy)
line(points3x,points3y,'Color','r')
hold off
grid on
axis equal
title("Move Pivot Points: O1-A-B-O2 colinear")

%% Define line equations
% Line 1 is the line containing points O11-A1-B1-O21 along a straight line.
Point1 = Pos_O1_Cayley;
dir1 = (Pos_P-Pos_A)/norm(Pos_P-Pos_A);

% Line 2 is the line containing points O12-A2-B2-O22 along a straight line.
Point2 = Pos_O2_Cayley;
dir2 = (Pos_P-Pos_B)/norm(Pos_P-Pos_B);

% Line 3 is the line containing points B-P-B1 along a straight line.
Point3 = Pos_P;
dir3 = (Pos_P-Pos_B)/norm(Pos_P-Pos_B);

% Line 4 is the line containing points A-P-B2 along a straight line.
Point4 = Pos_P;
dir4 = (Pos_P-Pos_A)/norm(Pos_P-Pos_A);

% Line 5 is the line containing points A1-P-A2 along a straight line.
Point5 = Pos_P;
dir5 = (Pos_A-Pos_B)/norm(Pos_A-Pos_B);

% Define known points of the Cayley diagram
Pos_O11_Cayley = Pos_O1_Cayley;
Pos_O12_Cayley = Pos_O2_Cayley;

% Find Intersection of Line L1 and L2
bvec = Point1(1:2) - Point2(1:2);
Amat = [-dir1(1:2) dir2(1:2)];
mres = Amat\bvec;
Pos_O21_Cayley = Point1 + mres(1)*dir1;

% Find Intersection of Line L1 and L3
bvec = Point1(1:2) - Point3(1:2);
Amat = [-dir1(1:2) dir3(1:2)];
mres = Amat\bvec;
Pos_B1_Cayley = Point1 + mres(1)*dir1;

% Find Intersection of Line L2 and L4
bvec = Point2(1:2) - Point4(1:2);
Amat = [-dir2(1:2) dir4(1:2)];
mres = Amat\bvec;
Pos_B2_Cayley = Point2 + mres(1)*dir2;

% Find Intersection of Line L1 and L5
bvec = Point1(1:2) - Point5(1:2);
Amat = [-dir1(1:2) dir5(1:2)];
mres = Amat\bvec;
Pos_A1_Cayley = Point1 + mres(1)*dir1;

% Find Intersection of Line L2 and L5
bvec = Point2(1:2) - Point5(1:2);
Amat = [-dir2(1:2) dir5(1:2)];
mres = Amat\bvec;
Pos_A2_Cayley = Point2 + mres(1)*dir2;

Pos_O22_Cayley = Pos_O21_Cayley;
Pos_P_Cayley = Pos_P;
Pos_B_Cayley = Pos_B;
Pos_A_Cayley = Pos_A;

points4x = [Pos_O11_Cayley(1,1) Pos_A1_Cayley(1,1) Pos_B1_Cayley(1,1) Pos_P(1,1) Pos_A1_Cayley(1,1) Pos_B1_Cayley(1,1)  Pos_O21_Cayley(1,1)];
points4y = [Pos_O11_Cayley(2,1) Pos_A1_Cayley(2,1) Pos_B1_Cayley(2,1) Pos_P(2,1) Pos_A1_Cayley(2,1) Pos_B1_Cayley(2,1)  Pos_O21_Cayley(2,1)];
points5x = [Pos_O12_Cayley(1,1) Pos_A2_Cayley(1,1) Pos_B2_Cayley(1,1) Pos_P(1,1) Pos_A2_Cayley(1,1) Pos_B2_Cayley(1,1)  Pos_O22_Cayley(1,1)];
points5y = [Pos_O12_Cayley(2,1) Pos_A2_Cayley(2,1) Pos_B2_Cayley(2,1) Pos_P(2,1) Pos_A2_Cayley(2,1) Pos_B2_Cayley(2,1)  Pos_O22_Cayley(2,1)];
%Figure 3
figno = figno + 1;
figure(figno)
plot(points3x,points3y,'ro')
hold on
plot(points4x,points4y,'ro')
plot(points5x,points5y,'ro')
line(pointsx,pointsy,'Color','k','LineStyle','--')
line(points3x,points3y,'Color','r','LineStyle','-')
line(points4x,points4y,'Color','b','LineStyle','-')
line(points5x,points5y,'Color','g','LineStyle','-')
hold off
grid on
title("Draw parallel lines get intersections")

%% The paralleograms
points1x_Cayley = [Pos_O12_Cayley(1,1) Pos_A2_Cayley(1,1) Pos_P_Cayley(1,1) Pos_B_Cayley(1,1) Pos_O12_Cayley(1,1)];
points1y_Cayley = [Pos_O12_Cayley(2,1) Pos_A2_Cayley(2,1) Pos_P_Cayley(2,1) Pos_B_Cayley(2,1) Pos_O12_Cayley(2,1)];

points2x_Cayley = [Pos_O11_Cayley(1,1) Pos_A1_Cayley(1,1) Pos_P_Cayley(1,1) Pos_A_Cayley(1,1) Pos_O11_Cayley(1,1)];
points2y_Cayley = [Pos_O11_Cayley(2,1) Pos_A1_Cayley(2,1) Pos_P_Cayley(2,1) Pos_A_Cayley(2,1) Pos_O11_Cayley(2,1)];

points3x_Cayley = [Pos_O21_Cayley(1,1) Pos_B1_Cayley(1,1) Pos_P_Cayley(1,1) Pos_B2_Cayley(1,1) Pos_O21_Cayley(1,1)];
points3y_Cayley = [Pos_O21_Cayley(2,1) Pos_B1_Cayley(2,1) Pos_P_Cayley(2,1) Pos_B2_Cayley(2,1) Pos_O21_Cayley(2,1)];

figno = figno + 1;
figure(figno)
plot(points1x_Cayley,points1y_Cayley,'k*')
hold on
plot(points2x_Cayley,points2y_Cayley,'ko')
plot(points3x_Cayley,points3y_Cayley,'ko')
line(points1x_Cayley,points1y_Cayley,'Color','r','LineStyle','--')
line(points2x_Cayley,points2y_Cayley,'Color','r','LineStyle','--')
line(points3x_Cayley,points3y_Cayley,'Color','r','LineStyle','--')
hold off
grid on
axis equal
%%
uvec1 = (Pos_A1_Cayley - Pos_P)/norm(Pos_A1_Cayley - Pos_P);
uvec2 = cross([0;0;1],uvec1);
uvec2 = uvec2/norm(uvec2);
Amat = [uvec1(1:2) uvec2(1:2)];
bvec = Pos_B1_Cayley(1:2) - Pos_P(1:2);
mres = Amat\bvec;

% New position
Pos_O11 = Pos_O1;
Pos_A1 = Pos_O11 + Pos_A1_Cayley - Pos_O11_Cayley;

uvec1 = (Pos_A1 - Pos_P)/norm(Pos_A1 - Pos_P);
uvec2 = cross([0;0;1],uvec1);
uvec2 = uvec2/norm(uvec2);
Pos_B1 = Pos_P + mres(1)*uvec1 + mres(2)*uvec2;

uvec1 = (Pos_A2_Cayley - Pos_P)/norm(Pos_A2_Cayley - Pos_P);
uvec2 = cross([0;0;1],uvec1);
uvec2 = uvec2/norm(uvec2);
Amat = [uvec1(1:2) uvec2(1:2)];
bvec = Pos_B2_Cayley(1:2) - Pos_P(1:2);
mres = Amat\bvec;

% New position
Pos_O12 = Pos_O2;
Pos_A2 = Pos_O12 + Pos_A2_Cayley - Pos_O12_Cayley;

uvec1 = (Pos_A2 - Pos_P)/norm(Pos_A2 - Pos_P);
uvec2 = cross([0;0;1],uvec1);
uvec2 = uvec2/norm(uvec2);
Pos_B2 = Pos_P + mres(1)*uvec1 + mres(2)*uvec2;

Pos_O21 = Pos_B1 + (Pos_B2 - Pos_P);
Pos_O22 = Pos_B2 + (Pos_B1 - Pos_P);

points1x = [Pos_O12(1,1) Pos_A2(1,1) Pos_P(1,1) Pos_B(1,1) Pos_O12(1,1)];
points1y = [Pos_O12(2,1) Pos_A2(2,1) Pos_P(2,1) Pos_B(2,1) Pos_O12(2,1)];

points2x = [Pos_O11(1,1) Pos_A1(1,1) Pos_P(1,1) Pos_A(1,1) Pos_O11(1,1)];
points2y = [Pos_O11(2,1) Pos_A1(2,1) Pos_P(2,1) Pos_A(2,1) Pos_O11(2,1)];

points3x = [Pos_O21(1,1) Pos_B1(1,1) Pos_P(1,1) Pos_B2(1,1) Pos_O21(1,1)];
points3y = [Pos_O21(2,1) Pos_B1(2,1) Pos_P(2,1) Pos_B2(2,1) Pos_O21(2,1)];
%
points4x = [Pos_O1(1,1) Pos_A(1,1) Pos_B(1,1) Pos_P(1,1) Pos_A(1,1) Pos_B(1,1) Pos_O2(1,1)];
points4y = [Pos_O1(2,1) Pos_A(2,1) Pos_B(2,1) Pos_P(2,1) Pos_A(2,1) Pos_B(2,1) Pos_O2(2,1)];

points5x = [Pos_O11(1,1) Pos_A1(1,1) Pos_B1(1,1) Pos_P(1,1) Pos_A1(1,1) Pos_B1(1,1) Pos_O21(1,1)];
points5y = [Pos_O11(2,1) Pos_A1(2,1) Pos_B1(2,1) Pos_P(2,1) Pos_A1(2,1) Pos_B1(2,1) Pos_O21(2,1)];

points6x = [Pos_O12(1,1) Pos_A2(1,1) Pos_B2(1,1) Pos_P(1,1) Pos_A2(1,1) Pos_B2(1,1) Pos_O22(1,1)];
points6y = [Pos_O12(2,1) Pos_A2(2,1) Pos_B2(2,1) Pos_P(2,1) Pos_A2(2,1) Pos_B2(2,1) Pos_O22(2,1)];

figno = figno + 1;
figure(figno)
plot(points4x,points4y,'ro')
hold on
plot(points5x,points5y,'ro')
plot(points6x,points6y,'ro')
line(points4x,points4y,'Color','k','LineStyle','-')
line(points5x,points5y,'Color','r','LineStyle','-')
line(points6x,points6y,'Color','b','LineStyle','-')
hold off
grid on
axis equal
title("The Cognates")
%% Output data
Pos_P1 = Pos_P;
Pos_P2 = Pos_P;
lg1 = norm(Pos_O11-Pos_O21);
len_O11_A1 = norm(Pos_O11-Pos_A1);
len_A1_B1 = norm(Pos_A1-Pos_B1);
len_O21_B1 = norm(Pos_O21-Pos_B1);
len_A1_P1 = norm(Pos_A1-Pos_P1);
len_B1_P1 = norm(Pos_B1-Pos_P1);

lg2 = norm(Pos_O12-Pos_O22);
len_O12_A2 = norm(Pos_O12-Pos_A2);
len_A2_B2 = norm(Pos_A2-Pos_B2);
len_O22_B2 = norm(Pos_O22-Pos_B2);
len_A2_P2 = norm(Pos_A2-Pos_P2);
len_B2_P2 = norm(Pos_B2-Pos_P2);

Pos_E1 = 0.5*(Pos_A1 + Pos_B1);
Pos_F1_N = Pos_P1 - Pos_E1;
dirx = (Pos_B1-Pos_A1)/norm(Pos_B1-Pos_A1);
diry = cross([0;0;1],dirx);
diry = diry/norm(diry);
Fij = [dirx(1:2) diry(1:2)]*Pos_F1_N(1:2);

[lg O1_A A_B O2_B A_P B_P]
[lg1 len_O11_A1 len_A1_B1 len_O21_B1 len_A1_P1 len_B1_P1]
[lg2 len_O12_A2 len_A2_B2 len_O22_B2 len_A2_P2 len_B2_P2]


