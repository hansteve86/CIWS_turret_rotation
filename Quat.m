%% Quaternion Project-Quaternion
clear ; close all
nx = [1, 0, 0]; %initial n-frame
ny = [0, 1, 0]; %initial n-frame
nz = [0, 0, 1]; %initial n-frame
n0 = [0, 0, 0]; %origin of n-frame

Theta1 = 0; %x this axis is our barrel
Theta2 = 90; %y
Theta3 = 180; %z

 % epsilon values
ep0 = cosd(Theta1/2) * cosd(Theta2/2) * cosd(Theta3/2) - sind(Theta1/2) * sind(Theta2/2) * sind(Theta3/2); 
ep1 = sind(Theta1/2) * cosd(Theta2/2) * cosd(Theta3/2) + sind(Theta2/2) * sind(Theta3/2) * cosd(Theta1/2);
ep2 = sind(Theta2/2) * cosd(Theta1/2) * cosd(Theta3/2) - sind(Theta1/2) * sind(Theta3/2) * cosd(Theta2/2);
ep3 = sind(Theta1/2) * sind(Theta2/2) * cosd(Theta3/2) + sind(Theta3/2) * cosd(Theta1/2) * cosd(Theta2/2);

 %epsilon vector consists of ep1 ep2 ep3
ep_vec = [ep1, ep2, ep3];
ThetaFinal = asind((ep1^2 + ep2^2 + ep3^2)^0.5)*2;
Lambda = ep_vec./sind(ThetaFinal/2); % finding unit vector Lambda

% finding new b frame
Bx = cosd(ThetaFinal)*nx + (1-cosd(ThetaFinal))*dot(nx,Lambda)*Lambda - sind(ThetaFinal)*(cross(nx, Lambda));
By = cosd(ThetaFinal)*ny + (1-cosd(ThetaFinal))*dot(ny,Lambda)*Lambda - sind(ThetaFinal)*(cross(ny, Lambda));
Bz = cosd(ThetaFinal)*nz + (1-cosd(ThetaFinal))*dot(nz,Lambda)*Lambda - sind(ThetaFinal)*(cross(nz, Lambda));

Steps = 120;
Thetas = linspace(0,ThetaFinal,Steps);
NFrame = [1 0 0; 0 1 0; 0 0 1];

TraceValues = zeros(2,Steps);
 %% Animation
FPS = 10;
LoopDelayFrames = 30;
Video = VideoWriter("Animation1");

Video.FrameRate = FPS;
open(Video);  
figure(1)

for Index = 1:Steps
    clf; %clear the plot
    
    Theta = Thetas(Index);
    Origin = num2cell([0 0 0]);
    [OriginX,OriginY,OriginZ] = deal(Origin{:});
    
    Nx = num2cell(NFrame(1,:));
    Ny = num2cell(NFrame(2,:));
    Nz = num2cell(NFrame(3,:));
    [Nxx,Nxy,Nxz] = deal(Nx{:});
    [Nyx,Nyy,Nyz] = deal(Ny{:});
    [Nzx,Nzy,Nzz] = deal(Nz{:});
    
    quiver3(OriginX,OriginY,OriginZ,Nxx,Nxy,Nxz,'r--', 'Linewidth', 4);
    text(1.1,0,0,"Nx","HorizontalAlignment","Center","FontSize",10);
    hold on
    quiver3(OriginX,OriginY,OriginZ,Nyx,Nyy,Nyz,'g--', 'Linewidth', 4);
    text(0,1.1,0,"Ny","HorizontalAlignment","Center","FontSize",10);
    hold on
    quiver3(OriginX,OriginY,OriginZ,Nzx,Nzy,Nzz,'b--', 'Linewidth', 4);
    text(0,0,1.1,"Nz","HorizontalAlignment","Center","FontSize",10);
    hold on
    
    CurrLambda = normalize(Lambda)*Theta * (1/100); %scale it down so we can see
    LambdaAxis = normalize(Lambda) * 10;
    plot3([-LambdaAxis(1),LambdaAxis(1)],[-LambdaAxis(2),LambdaAxis(2)],[-LambdaAxis(3),LambdaAxis(3)],'m--', 'Linewidth', 1);
    hold on
    quiver3(OriginX,OriginY,OriginZ,CurrLambda(1),CurrLambda(2),CurrLambda(3),'m', 'Linewidth', 2);
    text(CurrLambda(1)*1.1,CurrLambda(2)*1.1,CurrLambda(3)*1.1,sprintf("%0.1f°",Theta),"HorizontalAlignment","Center","FontSize",10);
    hold on
      
    nRb = [
        cosd(Theta)*nx + (1-cosd(Theta))*dot(nx,Lambda)*Lambda - sind(Theta)*(cross(nx, Lambda));
        cosd(Theta)*ny + (1-cosd(Theta))*dot(ny,Lambda)*Lambda - sind(Theta)*(cross(ny, Lambda));
        cosd(Theta)*nz + (1-cosd(Theta))*dot(nz,Lambda)*Lambda - sind(Theta)*(cross(nz, Lambda))
    ];
    BarrelX = nRb(1,1)*2;
    BarrelY = nRb(1,2)*2;
    BarrelZ = nRb(1,3)*2;
    quiver3(OriginX,OriginY,OriginZ,BarrelX,BarrelY,BarrelZ,'r', 'Linewidth', 4);
    text(BarrelX*1.1,BarrelY*1.1,BarrelZ*1.1,"Bx(Barrel)","HorizontalAlignment","Center","FontSize",16);
    hold on
    plot3( [0 0 nRb(2,1)], [0 0 nRb(2,2)], [0 0 nRb(2,3)], 'g', 'Linewidth', 2);
    text(nRb(2,1)*1.1,nRb(2,2)*1.1,nRb(2,3)*1.1,"By","HorizontalAlignment","Center","FontSize",10);
    hold on
    plot3( [0 0 nRb(3,1)], [0 0 nRb(3,2)], [0 0 nRb(3,3)], 'b', 'Linewidth', 2); 
    text(nRb(3,1)*1.1,nRb(3,2)*1.1,nRb(3,3)*1.1,"Bz","HorizontalAlignment","Center","FontSize",10);
    hold on
  
    TraceValues(:,Index) = [BarrelX BarrelY];

    title("Quaternions");
    axis([-2 2  -2 2  -2 2])
    %view(45,45);
    CurrentFigure = gcf; %reference to the current figure
    Frame = getframe(CurrentFigure);%"capture" the current figure as a frame
    writeVideo(Video,Frame); %write the frame to our video
end
for n = 1:LoopDelayFrames %keep writing the same figure for a few more seconds
    pause(1/FPS)
    CurrentFigure = gcf;
    Frame = getframe(CurrentFigure);
    writeVideo(Video,Frame);
end

hold off
figure(3)
Video2 = VideoWriter("Animation2");
Video2.FrameRate = FPS;
open(Video2);  
    
for Index = 1:Steps
    %clf; %clear the plot
    if Index > 1
        PrevTraceX = TraceValues(1,Index - 1);
        PrevTraceY = TraceValues(2,Index - 1);
        TraceX = TraceValues(1,Index);
        TraceY = TraceValues(2,Index);
        plot([PrevTraceX TraceX],[PrevTraceY TraceY],'c', 'Linewidth', 2)
        hold on
        
        title("Quaternions XY Trace");
        axis([-2 2 -2 2])
        axis square
        CurrentFigure = gcf; %reference to the current figure
        Frame = getframe(CurrentFigure);%"capture" the current figure as a frame
        writeVideo(Video2,Frame); %write the frame to our video   
    end
end

for n = 1:LoopDelayFrames %keep writing the same figure for a few more seconds
    pause(1/FPS)
    CurrentFigure = gcf;
    Frame = getframe(CurrentFigure);
    writeVideo(Video2,Frame);
end
