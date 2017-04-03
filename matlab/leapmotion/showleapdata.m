% showleapdata.m
%
% Show the fingertips 3d position feature of Leap Motion.

% Giulio Marin
%
% 2014/11/15
% giulio.marin@me.com

close all
clear

%% Parameters

filepath = 'leap_data.csv';
rowFingertips3d = 5;
rowHandDirection = 6;
rowPalmNormal = 9;
rowPalmPosition = 11;

%% Load data

leapData = csvread(filepath, 0, 1);

% Swap y and z axis
swap = [1 3 2];

palmPosition = leapData(rowPalmPosition,1:3)';
palmPosition = palmPosition(swap, :);
fingertips = reshape(leapData(rowFingertips3d,1:15),5,3)';
fingertips = fingertips(swap, :);
fingertips( :, ~any(fingertips,1) ) = []; % remove null columns
fingertips = fingertips - repmat(palmPosition, 1, size(fingertips,2));
palmNormal = leapData(rowPalmNormal,1:3)';
palmNormal = palmNormal(swap, :);
handDirection = leapData(rowHandDirection,1:3)';
handDirection = handDirection(swap, :);

palmPosition = palmPosition - palmPosition;

% move fingertips away from the plane
enhanceFactor = 0;
fingertips = fingertips + palmNormal * dot(fingertips - repmat(palmPosition, 1, size(fingertips,2)), repmat(palmNormal, 1, size(fingertips,2)));

%% Plane and projections

% Flip palmNormal
if palmNormal(3) > 0
    palmNormal = -palmNormal;
end

% Assign coefficients
a = palmNormal(1);
b = palmNormal(2);
c = palmNormal(3);
d = -[a b c] * palmPosition;

% Plane
nPoints = 2;
xAdd = 10;
yAdd = 50;
[xPlane,yPlane] = meshgrid(linspace(min(fingertips(1,:)) - xAdd, max(fingertips(1,:)) + xAdd,nPoints), linspace(min(fingertips(2,:)) - yAdd/2, max(fingertips(2,:)) + yAdd, nPoints));
zPlane = -(a*xPlane + b*yPlane + d) / c;

% Axis
l = 50;
n = palmPosition + l * palmNormal;
h = palmPosition + l * handDirection;

% Projected points
projectedFingertips = fingertips - palmNormal * dot(fingertips - repmat(palmPosition, 1, size(fingertips,2)), repmat(palmNormal, 1, size(fingertips,2)));

%% Show results

figure
hold on
grid on
axis equal

% Axis
hn = plot3([palmPosition(1) n(1)],[palmPosition(2) n(2)],[palmPosition(3) n(3)],'r','LineWidth',3);
hh = plot3([palmPosition(1) h(1)],[palmPosition(2) h(2)],[palmPosition(3) h(3)],'g','LineWidth',3);

% Plane
hp = surf(xPlane(1,:), yPlane(:,1), zPlane, 'FaceColor', [0.5 0.5 0.5]);
alpha(0.5);

% Fingertips and palm
hpc = scatter3(palmPosition(1,:),palmPosition(2,:),palmPosition(3,:),350,'fill','r');
hf = scatter3(fingertips(1,:),fingertips(2,:),fingertips(3,:),300,'fill','b');

% Projections
hproj = scatter3(projectedFingertips(1,:),projectedFingertips(2,:),projectedFingertips(3,:),50,'fill','k');

% Lines
hdist = plot3([palmPosition(1) fingertips(1,1)],[palmPosition(2) fingertips(2,1)],[palmPosition(3) fingertips(3,1)],'k','LineWidth',1.5);
helev = plot3([projectedFingertips(1,1) fingertips(1,1)],[projectedFingertips(2,1) fingertips(2,1)],[projectedFingertips(3,1) fingertips(3,1)],'--k','LineWidth',1.5);
for i=1:size(fingertips,2)
    plot3([palmPosition(1) fingertips(1,i)],[palmPosition(2) fingertips(2,i)],[palmPosition(3) fingertips(3,i)],'k','LineWidth',1.5)
    plot3([palmPosition(1) projectedFingertips(1,i)],[palmPosition(2) projectedFingertips(2,i)],[palmPosition(3) projectedFingertips(3,i)],':k','LineWidth',1.5)
    plot3([projectedFingertips(1,i) fingertips(1,i)],[projectedFingertips(2,i) fingertips(2,i)],[projectedFingertips(3,i) fingertips(3,i)],'--k','LineWidth',1.5)
end

% Set visualization
xlabel('x','FontSize',20);ylabel('z','FontSize',20);zlabel('y','FontSize',20)

set(gca,'yDir','reverse')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])

legend([hn,hh,hp,hpc,hf,hproj,hdist,helev],{'Palm normal n', 'Hand direction h', 'Palm plane', 'Palm center', 'Fingertip 3D positions', 'Fingertip projections', 'Fingertip distances', 'Fingertip elevations'})
hold off