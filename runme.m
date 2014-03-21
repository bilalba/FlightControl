clear all
close all
clc

figure;
axis equal;
axis([0 100 0 100]);
hold on;
ballPos = [50 50];
ballRad = 4;
pathPos(1,:)=ballPos;
x=plot(pathPos(:,1),pathPos(:,2), 'linewidth', 3);
r = rectangle('position', [ballPos-ballRad 2*ballRad*[1 1]],...
    'curvature', [1 1], 'facecolor', 'r');
e4=tic;
mouseLoc = zeros(2,3);

pressed = 0;
justClicked = 0;
set(gcf, 'WindowButtonMotionFcn', 'mouseLoc = get(gca, ''CurrentPoint'');',...
    'WindowButtonDownFcn', 'pressed = 1; justClicked = 1;',...
    'WindowButtonUpFcn', 'pressed = 0;');

xyz=0
dr=0;
sn = 1;
dragging = 0;
i=1;
aw=0
bb=1



while 1
    
    start=tic;
    
    if pressed
        
        if justClicked
            justClicked = 0;
            pathPos=[0 0]
            pathPos(1,:)=ballPos;
            i=1;
            if norm(mouseLoc(1,1:2)-ballPos) <= ballRad
                dragging = 1;
            end
        end
        if dragging
            
            dr=1;
           
            
            
            if norm(mouseLoc(1,1:2)-pathPos(i,:))~=0
                i=i+1;
                
                pathPos(i,:) = mouseLoc(1,1:2);
                
                
                
                set(x, 'Xdata',pathPos(:,1), 'Ydata', pathPos(:,2))
                
                
                
            end
            
            
            
        end
    else
        dragging = 0;
    end
    if dr==1 && (i)~=1
        direct=[pathPos(2,1)-pathPos((1),1) (pathPos(2,2)-pathPos((1),2))];
        grad=direct/norm(direct);
        
        
        
        xCord=grad(1);
        yCord=grad(2);
        
        
        distance=norm(pathPos(2,1:2)-pathPos(1,1:2));

        veloC=5;
        
        
      
        
        veloC=veloC*0.01
        
         iter=floor(distance/veloC)
         rem=distance-veloC*iter
         extra=rem/iter
         
         
         veloC=veloC+extra
         
           time=(distance/veloC)*0.01
        
        xCord=xCord*veloC
        yCord=yCord*veloC
        ballPos = [ballPos(1)+xCord ballPos(2)+yCord];
        
        
        
        set(r, 'position', [ballPos-ballRad 2*ballRad*[1 1]])
        aw=aw+1
        x123=toc(e4)
        x123=x123-toc(start)
        if iter<=aw
            pathPos=pathPos(2:(i),:)
            i=i-1;
            pathPos
            
            e4=tic
            
            set(x, 'Xdata',pathPos(:,1), 'Ydata', pathPos(:,2))
            aw=0
            
            xyz=1
        end
    elseif xyz==1
        
        ballPos = [ballPos(1)+xCord ballPos(2)+yCord];
        
        pathPos(1,:)=ballPos;
        
        set(r, 'position', [ballPos-ballRad 2*ballRad*[1 1]])
        
        
    end
    toc(start)
    
    if xyz==1
    title(180*(atan2(yCord,xCord))/pi+(180))
    end
    pause(0.01-toc(start));
    
end