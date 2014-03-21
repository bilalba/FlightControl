clear all
close all
clc

figure;
axis equal;
axis([0 100 0 100]);
hold on;
ballPos = [50 50];
ballRad = 4;
r = rectangle('position', [ballPos-ballRad 2*ballRad*[1 1]],...
    'curvature', [1 1], 'facecolor', 'r');
e4=tic;
mouseLoc = zeros(2,3);

pressed = 0;
justClicked = 0;
set(gcf, 'WindowButtonMotionFcn', 'mouseLoc = get(gca, ''CurrentPoint'');',...
    'WindowButtonDownFcn', 'pressed = 1; justClicked = 1;',...
    'WindowButtonUpFcn', 'pressed = 0;');
dr=0;
sn = 1;
dragging = 0;
i=1;
u=2;
aw=0
bb=1
pathPos(1,:,bb)=ballPos;

while 1
    
    start=tic;
    
    if pressed
        
        if justClicked
            justClicked = 0;
            if norm(mouseLoc(1,1:2)-ballPos) <= ballRad
                dragging = 1;
            end
        end
        if dragging
            
            dr=1;
           
            
            
            if norm(mouseLoc(1,1:2)-pathPos(i,:))~=0
                i=i+1;
                
                pathPos(i,:) = mouseLoc(1,1:2);
                
                eval(['plot(pathPos(:,1),pathPos(:,2), 'linewidth', 3)' ])
                
                
                
                
            end
            
            
            
        end
    else
        dragging = 0;
    end
    if dr==1 && (u-i)~=0
        direct=[pathPos(u,1)-pathPos((u-1),1) (pathPos(u,2)-pathPos((u-1),2))];
        grad=direct/norm(direct);
        
        
        
        xCord=grad(1);
        yCord=grad(2);
        
        
        distance=norm(pathPos(u,1:2)-pathPos(u-1,1:2))

        veloC=5;
        
        
      
        
        veloC=veloC*0.04
        
         iter=floor(distance/veloC)
         rem=distance-veloC*iter
         extra=rem/iter
         
         
         veloC=veloC+extra
         
           time=(distance/veloC)*0.04
        
        xCord=xCord*veloC
        yCord=yCord*veloC
        ballPos = [ballPos(1)+xCord ballPos(2)+yCord];        
        set(r, 'position', [ballPos-ballRad 2*ballRad*[1 1]])
        aw=aw+1
        x123=toc(e4)
        x123=x123-toc(start)
        if iter<=aw
            
            u=u+1
            e4=tic
            aw=0
        end
    end
    
   
    pause(0.04-toc(start));
    
end