%% Quaternion Project-Euler
clear ; close all

Steps = 120;
ZAngle = 180;
YAngle = -90;

SplitReferenceFrames = 0; % 1 = true, 0 = false;
FrameRatio = 0.5;%ratio pf steps to give each reference frame A;

ZSteps = Steps;
YSteps = Steps;

if SplitReferenceFrames == 1 
    ZSteps = floor(Steps*FrameRatio);
    YSteps = ceil(Steps* (1-FrameRatio));
end

ZRot = linspace(0,ZAngle,ZSteps); %Yaw
YRot = linspace(0,YAngle,YSteps); %Pitch
NFrame = [1 0 0; 0 1 0; 0 0 1];

TraceValues = zeros(2,Steps);
 %% Animation
FPS = 10;
LoopDelayFrames = 30;
Video = VideoWriter("Animation1");
Video.FrameRate = FPS;
open(Video);  
figure(1);

for Index = 1:Steps
    clf; %clear the plot
    if SplitReferenceFrames == 1 
        Theta1 = ZAngle;
        if Index <= ZSteps
            Theta1 = ZRot(Index);
        end
        Theta2 = 0;
        if Index > ZSteps
            Theta2 = YRot(Index - ZSteps);
        end
    else
        Theta1 = ZRot(Index);
        Theta2 = YRot(Index);
    end
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
    hold on
    
    nRa = [cosd(Theta1), sind(Theta1), 0; -sind(Theta1), cosd(Theta1), 0; 0, 0, 1];
    
    plot3([0 0 nRa(1,1)], [0 0 nRa(1,2)], [0 0 nRa(1,3)], 'r', 'Linewidth', 2);
    text(nRa(1,1)*1.1,nRa(1,2)*1.1,nRa(1,3)*1.1,"Ax","HorizontalAlignment","Center","FontSize",10);
    hold on
    plot3( [0 0 nRa(2,1)], [0 0 nRa(2,2)], [0 0 nRa(2,3)], 'g', 'Linewidth', 2);
    text(nRa(2,1)*1.1,nRa(2,2)*1.1,nRa(2,3)*1.1,"By = Ay","HorizontalAlignment","Center","FontSize",10);
    hold on
    plot3( [0 0 nRa(3,1)], [0 0 nRa(3,2)], [0 0 nRa(3,3)], 'b', 'Linewidth', 2);
    text(nRa(3,1)*1.1,nRa(3,2)*1.1,nRa(3,3)*1.1,"Nz = Az","HorizontalAlignment","Center","FontSize",10);
    grid on
    
    
    aRb = [cosd(Theta2), 0, -sind(Theta2); 0, 1, 0; sind(Theta2), 0, cosd(Theta2)];
    nRb = aRb*nRa;
    BarrelX = nRb(1,1)*2;
    BarrelY = nRb(1,2)*2;
    BarrelZ = nRb(1,3)*2;
    quiver3(OriginX,OriginY,OriginZ,BarrelX,BarrelY,BarrelZ,'r', 'Linewidth', 4);
    text(BarrelX*1.1,BarrelY*1.1,BarrelZ*1.1,"Bx(Barrel)","HorizontalAlignment","Center","FontSize",16);
    hold on
    plot3( [0 0 nRb(2,1)], [0 0 nRb(2,2)], [0 0 nRb(2,3)], 'g', 'Linewidth', 2);
    hold on
    plot3( [0 0 nRb(3,1)], [0 0 nRb(3,2)], [0 0 nRb(3,3)], 'b', 'Linewidth', 2); 
    text(nRb(3,1)*1.1,nRb(3,2)*1.1,nRb(3,3)*1.1,"Bz","HorizontalAlignment","Center","FontSize",10);
    hold on
    
    TraceValues(:,Index) = [BarrelX BarrelY];
    
    if SplitReferenceFrames == 1 
        title("Euler Angles nRa + aRb");
    else
        title("Euler Angles nRb");
    end
    
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
        
        if SplitReferenceFrames == 1 
            title("Euler Angles nRa + aRb XY Trace");
        else
            title("Euler Angles nRb XY Trace");
        end
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