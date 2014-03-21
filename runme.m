%% Initilization

clear all
close all
clc

figure;
axis equal;
axis([0 200 0 100]);
hold on;


im = imread('screen.jpg');  % read image; 
im = imrotate(im,360)

imagesc([0 200], [0 100], flipdim(im,1));

set(gca,'Ydir','normal');

%% Variables
ts=1
sn = 1;
ballRad = 2;
e4=tic;
num=2;
for bb=1:num
    distance(bb)=0;
     veloC(bb)=5;
     veloC(bb)=veloC(bb)*0.02;
    iter(bb)=0;
    dragging(bb) = 0;
    dr(bb)=0;
    i(bb)=1;
    aw(bb)=0;
    xyz(bb)=0;
    
    %% BALL POS
    
    if rand<0.25
        ballPos(bb,:) = [round(200*rand) 4];
        angel=rand*pi
        xCord(bb)=cos(angel)*veloC(bb);
        yCord(bb)=sin(angel)*veloC(bb);
        
    elseif rand<0.5
        ballPos(bb,:) = [round(200*rand) 96];
        angel=-rand*pi
        xCord(bb)=cos(angel)*veloC(bb);
        yCord(bb)=sin(angel)*veloC(bb);
        
    elseif rand<0.75
    ballPos(bb,:) = [4 round(100*rand)];
    angel=-rand*pi+pi/2;
        xCord(bb)=cos(angel)*veloC(bb);
        yCord(bb)=sin(angel)*veloC(bb);
      
    else
        ballPos(bb,:) = [196 round(100*rand)];
        angel=-rand*pi+pi/2;
        xCord(bb)=cos(angel)*veloC(bb);
        yCord(bb)=sin(angel)*veloC(bb);
    end
    
    %%
    eval(['pathPos_' num2str(bb) '(1,1:2)=ballPos(bb,:)']);
    
    %% PATH
    
    eval(['plot_' num2str(bb) '=plot(pathPos_' num2str(bb) '(:,1),pathPos_' num2str(bb) '(:,2)), ''linewidth'', 5']);
                    
    
    %% Ball
  
    r(bb) = rectangle('position', [ballPos(bb,:)-ballRad 2*ballRad*[1 1]],...
        'curvature', [1 1], 'facecolor', 'r');
    
    
    
    
    
    
end

%% Mouse code

mouseLoc = zeros(2,3);

pressed = 0;
justClicked = 0;
set(gcf, 'WindowButtonMotionFcn', 'mouseLoc = get(gca, ''CurrentPoint'');',...
    'WindowButtonDownFcn', 'pressed = 1; justClicked = ones(num);',...
    'WindowButtonUpFcn', 'pressed = 0;');

    %% Start of main loop



while ts
    
    start=tic;
    for bb=1:num
        
        if pressed
            
            
            if justClicked(bb)
                justClicked(bb) = 0;
                
                
                
                
                if norm(mouseLoc(1,1:2)-ballPos(bb,:)) <= ballRad
                    
                    eval(['pathPos_' num2str(bb) '=[0 0]']);
                eval(['pathPos_' num2str(bb) '(1,:)=ballPos(bb,:)']);
                i(bb)=1;
                    
                    dragging(bb) = 1
                    
                end
                
            end
            
            
            if dragging(bb)
                
                
                
                
                dr(bb)=1;
                
                
                frt=mouseLoc(1,1:2);
                
                  for rt=1:2
                       
                        if  frt(1)>195
                            frt(1)=195;
                        end
                        if frt(2)>95
                           frt(2)=95; 
                        end
                        
                        if frt(rt)<5
                           frt(rt)=5;
                        end
                    end
                
                
                if eval(['norm(frt-pathPos_' num2str(bb) '(i(bb),1:2))~=0'])
                    i(bb)=i(bb)+1
                    
                    eval(['pathPos_' num2str(bb) '(i(bb),:)=frt']);
                    
                  
                    
                    
                    eval(['set(plot_' num2str(bb) ', ''Xdata'', pathPos_' num2str(bb) '(:,1), ''Ydata'', pathPos_' num2str(bb) '(:,2))']);
                    
                              
                    
                    
                    
                end
                
            end
            
            
            
        else
            dragging(bb) = 0;
        end
    end
    for bb=1:num
        
        if dr(bb)==1 && i(bb)~=1
            eval(['direct(bb,:)=[pathPos_' num2str(bb) '(2,1)-pathPos_' num2str(bb) '(1,1) pathPos_' num2str(bb) '(2,2)-pathPos_' num2str(bb) '(1,2)]'])
            
            
            grad(bb,:)=direct(bb,:)/norm(direct(bb,:))
            
            
            
            xCord(bb)=grad(bb,1)
            yCord(bb)=grad(bb,2)
            
            eval(['distance(bb)=norm(pathPos_' num2str(bb) '(2,1:2)-pathPos_' num2str(bb) '(1,1:2))'])
            
           
            
            veloC(bb)=5;
     veloC(bb)=veloC(bb)*0.02;
            
           
            
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
                
                
                eval(['pathPos_' num2str(bb) '=pathPos_' num2str(bb) '(2:(i(bb)),:)']);
                
                i(bb)=i(bb)-1;
                
                 eval(['set(plot_' num2str(bb) ', ''Xdata'', pathPos_' num2str(bb) '(:,1), ''Ydata'', pathPos_' num2str(bb) '(:,2))']);
                    
                aw(bb)=0
                
                xyz(bb)=1
            end
            
        elseif xyz(bb)==1
            
            ballPos(bb,:)=[ballPos(bb,1)+xCord(bb) ballPos(bb,2)+yCord(bb)];
            
            eval(['pathPos_' num2str(bb) '(1,:)=ballPos(bb,:)']);
            
            set(r(bb), 'position', [ballPos(bb,:)-ballRad 2*ballRad*[1 1]])
            
        else
             ballPos(bb,:)=[ballPos(bb,1)+xCord(bb) ballPos(bb,2)+yCord(bb)];
             
             set(r(bb), 'position', [ballPos(bb,:)-ballRad 2*ballRad*[1 1]])
            
            
        end
    
    %% Collision    
        
    uv=1:num;
    uv(bb)=[];
    
    
    for uv=uv
        
       if norm(ballPos(bb,:)-ballPos(uv,:))<4
          
           ts=0;
       end
        
    end
    
    
    if ballPos(bb,1)>196 || ballPos(bb,1)<4
        xCord(bb)=-xCord(bb)
    end
    
    if ballPos(bb,2)>96 || ballPos(bb,2)<4
        yCord(bb)=-yCord(bb)
    end
        
    %% LANDING
    
    if eval(['norm(pathPos_' num2str(bb) '(end,:)-[50,50])<5'])
    
        dragging(bb)=0;
        
    end
    end
    
    if rand>0.999
        
        num=num+1;
        
        veloC(num)=5;
     veloC(num)=veloC(num)*0.02;
    
    distance(num)=0;
 
    iter(num)=0;
    dragging(num) = 0;
    dr(num)=0;
    i(num)=1;
    aw(num)=0;
    xyz(num)=0;
    
    %% BALL POS
    
    if rand<0.25
        ballPos(num,:) = [round(200*rand) 4];
        angel=rand*pi
        xCord(num)=cos(angel)*veloC(num);
        yCord(num)=sin(angel)*veloC(num);
        
    elseif rand<0.5
        ballPos(num,:) = [round(200*rand) 96];
        angel=-rand*pi
        xCord(num)=cos(angel)*veloC(num);
        yCord(num)=sin(angel)*veloC(num);
        
    elseif rand<0.75
    ballPos(num,:) = [4 round(100*rand)];
    angel=-rand*pi+pi/2;
        xCord(num)=cos(angel)*veloC(num);
        yCord(num)=sin(angel)*veloC(num);
      
    else
        ballPos(num,:) = [196 round(100*rand)];
        angel=-rand*pi+pi/2;
        xCord(num)=cos(angel)*veloC(num);
        yCord(num)=sin(angel)*veloC(num);
    end
    
    %%
    eval(['pathPos_' num2str(num) '(1,1:2)=ballPos(num,:)']);
    
    %% PATH
    
    eval(['plot_' num2str(num) '=plot(pathPos_' num2str(num) '(:,1),pathPos_' num2str(num) '(:,2)), ''linewidth'', 5']);
                    
    
    %% Ball
  
    r(num) = rectangle('position', [ballPos(num,:)-ballRad 2*ballRad*[1 1]],...
        'curvature', [1 1], 'facecolor', 'r');
    
    
    end
    

    
    
    
    
   pause(0.02-toc(start));
    
end