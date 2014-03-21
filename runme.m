%% Initilization

clear all
close all
clc

figure;
axis equal;
axis([0 100 0 100]);
hold on;

%%

ballRad = 4;

e4=tic;
%% Mouse code

mouseLoc = zeros(2,3);

pressed = 0;
justClicked = 0;
set(gcf, 'WindowButtonMotionFcn', 'mouseLoc = get(gca, ''CurrentPoint'');',...
    'WindowButtonDownFcn', 'pressed = 1;',...
    'WindowButtonUpFcn', 'pressed = 0;');


%% Variables
sn = 1;
i=1;
u=2;
aw=0

for bb=1:2
    distance(bb)=0;
    justClicked(bb)=1;
    iter(bb)=0;
    dragging(bb) = 0;
    dr(bb)=0;
    i(bb)=1;
    u(bb)=2;
    aw(bb)=0
    
    %% Ball
    ballPos(bb,:) = round(100*rand(1,2))
    r(bb) = rectangle('position', [ballPos(bb,:)-ballRad 2*ballRad*[1 1]],...
        'curvature', [1 1], 'facecolor', 'r');
    eval(['pathPos_' num2str(bb) '(1,1:2)=ballPos(bb,:)']);
    
    
    
    
    
end

%% Start of main loop



while 1
    
    start=tic;
    
    for bb=1:2
        
        if pressed
            
            
            if justClicked(bb)
                justClicked(bb) = 0;
                
                bb
                if norm(mouseLoc(1,1:2)-ballPos(bb,:)) <= ballRad
                    dragging(bb) = 1
                    
                end
                
            end
            
            
            if dragging(bb)
                
                
                
                
                dr(bb)=1;
                
                if eval(['norm(mouseLoc(1,1:2)-pathPos_' num2str(bb) '(i(bb),1:2))~=0'])
                    i(bb)=i(bb)+1
                    
                    eval(['pathPos_' num2str(bb) '(i(bb),:)=mouseLoc(1,1:2)']);
                    
                    eval(['plot_' num2str(bb) '=plot(pathPos_' num2str(bb) '(:,1),pathPos_' num2str(bb) '(:,2))']);
                    
                    
                    
                    
                    
                end
                
            end
            
            
            
        else
            dragging(bb) = 0;
        end
    end
    for bb=1:2
        
        if dr(bb)==1 && (u(bb)-i(bb))~=0
            eval(['direct(bb,:)=[pathPos_' num2str(bb) '(u(bb),1)-pathPos_' num2str(bb) '((u(bb)-1),1) pathPos_' num2str(bb) '(u(bb),2)-pathPos_' num2str(bb) '((u(bb)-1),2)]'])
            
            
            grad(bb,:)=direct(bb,:)/norm(direct(bb,:))
            
            
            
            xCord(bb)=grad(bb,1)
            yCord(bb)=grad(bb,2)
            
            eval(['distance(bb)=norm(pathPos_' num2str(bb) '(u(bb),1:2)-pathPos_' num2str(bb) '(u(bb)-1,1:2))'])
            
           
            
            
            
            veloC(bb)=5;
            
            
            
            
            veloC(bb)=veloC(bb)*0.04;
            
            iter(bb)=floor(distance(bb)/veloC(bb))
            rem(bb)=distance(bb)-veloC(bb)*iter(bb)
            extra(bb)=rem(bb)/iter(bb)
            
            
            veloC(bb)=veloC(bb)+extra(bb)
            
            
            
            xCord(bb)=xCord(bb)*veloC(bb);
            yCord(bb)=yCord(bb)*veloC(bb);
            ballPos(bb,:) = [ballPos(bb,1)+xCord(bb) ballPos(bb,2)+yCord(bb)]
            set(r(bb), 'position', [ballPos(bb,:)-ballRad 2*ballRad*[1 1]])
            aw(bb)=aw(bb)+1
            
            if iter(bb)<=aw(bb)
                
                u(bb)=u(bb)+1
                
                aw(bb)=0
            end
        end
        
    end
    
    
    pause(0.04-toc(start));
    
end