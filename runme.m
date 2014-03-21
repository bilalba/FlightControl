%% Initilization




clear all
close all
clc

figure;
axis equal;
axis([0 200 0 100]);
hold on;
ruv=0.995

%% Start screen

gui = imread('gui.png');
gui1=imagesc([0 200], [0 100], flipdim(gui,1));

pause

%% Read images and audio
im = imread('screen.jpg');  % read image;
[alert,fs]= wavread('this.wav');
alert1=audioplayer(alert,fs);

[crash,fs1]= wavread('this1.wav');
crash1=audioplayer(crash,fs1);

[sprite map alphasprite]=imread('low_qual.png'); %read image sprite


[h1 w d] = size(sprite);

imagesc([0 200], [0 100], flipdim(im,1));

set(gca,'Ydir','normal');

%% Variables
ts=1;
sn = 1;
ballRad = 4;
e4=tic;
num=2;
zxca=0;
score=0;
ney=text(130,90,['Your score is:' num2str(score)], 'fontsize', 30, 'fontname', 'georgia');
   
distance=zeros(num,1);
     veloC=zeros(num,1);
    iter=zeros(num,1);
    dragging= zeros(num,1);
    dr=zeros(num,1);
    i=zeros(num,1);
    aw=zeros(num,1);
    xyz=zeros(num,1);
    landing=zeros(num,1);
    sv=zeros(num,1);
    zv=zeros(num,1);
    ballPos=zeros(num,2);
    

for bb=1:num
    distance(bb)=0;
     veloC(bb)=8;
     veloC(bb)=veloC(bb)*0.03;
    iter(bb)=0;
    dragging(bb) = 0;
    dr(bb)=0;
    i(bb)=1;
    aw(bb)=0;
    xyz(bb)=0;
    landing(bb)=0;
    sv(bb)=0;
    zv(bb)=0;
    
    %% Ball Position define
    
    if rand<0.25
        ballPos(bb,:) = [round(200*rand) 4];
        angel(bb)=rand*pi;
        xCord(bb)=cos(angel(bb))*veloC(bb);
        yCord(bb)=sin(angel(bb))*veloC(bb);
        
    elseif rand<0.5
        ballPos(bb,:) = [round(200*rand) 96];
        angel(bb)=-rand*pi;
        xCord(bb)=cos(angel(bb))*veloC(bb);
        yCord(bb)=sin(angel(bb))*veloC(bb);
        
    elseif rand<0.75
    ballPos(bb,:) = [4 round(100*rand)];
    angel(bb)=-rand*pi+pi/2;
        xCord(bb)=cos(angel(bb))*veloC(bb);
        yCord(bb)=sin(angel(bb))*veloC(bb);
      
    else
        ballPos(bb,:) = [196 round(100*rand)];
        angel(bb)=-rand*pi+pi/2;
        xCord(bb)=cos(angel(bb))*veloC(bb);
        yCord(bb)=sin(angel(bb))*veloC(bb);
    end
    
    [sv(bb) zv(bb)]=RoundOffAngle(angel(bb));
    
    
    
    
    %%
    eval(['pathPos_' num2str(bb) '(1,1:2)=ballPos(bb,:)']);
    
    %% PATH
    
    eval(['plot_' num2str(bb) '=plot(pathPos_' num2str(bb) '(:,1),pathPos_' num2str(bb) '(:,2), ''linewidth'', 2, ''lineStyle'', ''--'', ''Color'',[0.4 0.4 0.4])']);  %%
                    
    
    %% Aeroplane
 
    h(bb)=imagesc([ballPos(bb,1)-ballRad ballPos(bb,1)+ballRad], [ballPos(bb,2)-ballRad ballPos(bb,2)+ballRad], flipdim(sprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:),3), 'alphadata', alphasprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:));
      
    
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
                
                
                
                
                if norm(mouseLoc(1,1:2)-ballPos(bb,:)) <= ballRad  %Drag
                    
                    eval(['pathPos_' num2str(bb) '=[0 0]']);
                eval(['pathPos_' num2str(bb) '(1,:)=ballPos(bb,:)']);
                i(bb)=1;
                    landing(bb)= 0;
                    dragging(bb) = 1
                    
                end
                
            end
            
            
            if dragging(bb)
                
                
                
                
                dr(bb)=1;
                
                
                frt=mouseLoc(1,1:2);
                
                  
                       
                        if  frt(1)>195
                            frt(1)=195;
                        end
                        if frt(2)>95
                           frt(2)=95; 
                        end
                        for rt=1:2
                            if frt(rt)<5
                                frt(rt)=5;
                            end
                        end
                
                
                if eval(['norm(frt-pathPos_' num2str(bb) '(i(bb),1:2))~=0'])
                    i(bb)=i(bb)+1
                    
                    eval(['pathPos_' num2str(bb) '(i(bb),:)=frt']);
                    
                  
                    
                    
                    eval(['set(plot_' num2str(bb) ', ''Xdata'', pathPos_' num2str(bb) '(:,1), ''Ydata'', pathPos_' num2str(bb) '(:,2), ''lineStyle'', ''--'', ''Color'', [0.4 0.4 0.4])']);
                    
                              
                    
                    
                    
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
            
            
            
            xCord(bb)=grad(bb,1);
            yCord(bb)=grad(bb,2);
            
            
    [sv(bb) zv(bb)]=RoundOffAngle(atan2(yCord(bb),xCord(bb))); 
            
            eval(['distance(bb)=norm(pathPos_' num2str(bb) '(2,1:2)-pathPos_' num2str(bb) '(1,1:2))'])
            
           
            
            veloC(bb)=8;
     veloC(bb)=veloC(bb)*0.03;
            
           
            
            iter(bb)=floor(distance(bb)/veloC(bb))
            rem(bb)=distance(bb)-veloC(bb)*iter(bb)
            extra(bb)=rem(bb)/(iter(bb)+1)
            
            
            veloC(bb)=veloC(bb)+extra(bb)
            
          
            
            xCord(bb)=xCord(bb)*veloC(bb);
            yCord(bb)=yCord(bb)*veloC(bb);
            ballPos(bb,:) = [ballPos(bb,1)+xCord(bb) ballPos(bb,2)+yCord(bb)]
            
            set(h(bb), 'Xdata', [ballPos(bb,1)-ballRad ballPos(bb,1)+ballRad], 'Ydata', [ballPos(bb,2)-ballRad ballPos(bb,2)+ballRad], 'CData', flipdim(sprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:),3), 'alphadata', alphasprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:));
    
            aw(bb)=aw(bb)+1
            
            if iter(bb)<=aw(bb)
                
                
                eval(['pathPos_' num2str(bb) '=pathPos_' num2str(bb) '(2:(i(bb)),:)']);
                
                i(bb)=i(bb)-1;
                
                 eval(['set(plot_' num2str(bb) ', ''Xdata'', pathPos_' num2str(bb) '(:,1), ''Ydata'', pathPos_' num2str(bb) '(:,2))']);
                    
                aw(bb)=0;
                
                xyz(bb)=1;
            end
            
        elseif xyz(bb)==1
            
            ballPos(bb,:)=[ballPos(bb,1)+xCord(bb) ballPos(bb,2)+yCord(bb)];
            
            eval(['pathPos_' num2str(bb) '(1,:)=ballPos(bb,:)']);
            
            
           set(h(bb), 'Xdata', [ballPos(bb,1)-ballRad ballPos(bb,1)+ballRad], 'Ydata', [ballPos(bb,2)-ballRad ballPos(bb,2)+ballRad], 'CData', flipdim(sprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:),3), 'alphadata', alphasprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:));
    
        else
             ballPos(bb,:)=[ballPos(bb,1)+xCord(bb) ballPos(bb,2)+yCord(bb)];
             
            
            set(h(bb), 'Xdata', [ballPos(bb,1)-ballRad ballPos(bb,1)+ballRad], 'Ydata', [ballPos(bb,2)-ballRad ballPos(bb,2)+ballRad], 'CData', flipdim(sprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:),3), 'alphadata', alphasprite(sv(bb)+1:sv(bb)+100,zv(bb)+1:zv(bb)+100,:));
    
        end
    
    %% Collision    
        
    uv=1:num;
   
    uv(bb)=[];
    
    for uv=uv
        poi=norm(ballPos(bb,:)-ballPos(uv,:));
        
        if poi<12    % if distance less than 12
          play(alert1);
          
       
           
        end
       if poi<8 % if distance less than 8
          stop(alert1);
          play(crash1);
         
           ts=0;
       end
        
    end
    
    
    if ballPos(bb,1)>196 || ballPos(bb,1)<4
        xCord(bb)=-xCord(bb)
    end
    
    if ballPos(bb,2)>96 || ballPos(bb,2)<4
        yCord(bb)=-yCord(bb)
    end
    
    
    [sv(bb) zv(bb)]=RoundOffAngle(atan2(yCord(bb),xCord(bb)));
        
    %% IF LANDING
    
    if eval(['norm(pathPos_' num2str(bb) '(end,:)-[80,50])<2']) && landing(bb)==0 && dragging(bb)==1
        
       lmn=atan2(eval(['pathPos_' num2str(bb) '(end,2)-pathPos_' num2str(bb) '(end-3,2)']),eval(['pathPos_' num2str(bb) '(end,1)-pathPos_' num2str(bb) '(end-3,1)']))
    
        
        if lmn<0.785 && lmn>-0.785
            
        dragging(bb)=0;
        landing(bb)=1;
        end
   
    end
    
    
    if landing(bb)==1;
       %Landing goes here
        eval(['set(plot_' num2str(bb) ', ''Xdata'', pathPos_' num2str(bb) '(:,1), ''Ydata'', pathPos_' num2str(bb) '(:,2), ''Color'', [0 1 0] )']);
                    
        if norm(ballPos(bb,:)-[80,50])<3
        
            
           zxca=bb;
        
        end
        
    end
    
    
    end
    
    if zxca~=0
        
        score=score+1;
        
        
        delete(ney);

    ney=text(130,90,['Your score is:' num2str(score)], 'fontsize', 30, 'fontname', 'georgia');

         distance(zxca)=[]
        iter(zxca)=[];
    dragging(zxca) = [];
    dr(zxca)=[];
    i(zxca)=[];
    aw(zxca)=[];
    xyz(zxca)=[];
    landing(zxca)=[];
    sv(zxca)=[];
    zv(zxca)=[];
    ballPos(zxca,:)=[];
    angel(zxca)=[];
    xCord(zxca)=[];
    yCord(zxca)=[];
    
    
    
    delete(h(zxca))
    h(zxca)=[];
    
    if zxca==num
        eval(['pathPos_' num2str(zxca) '=[]'])
    else
    
    for gtr=zxca+1:num
    eval(['pathPos_' num2str(gtr-1) '=[]'])
    eval(['pathPos_' num2str(gtr-1) '=pathPos_' num2str(gtr)])
    
    end
    
    eval(['pathPos_' num2str(num) '=[]'])
    end
    num=num-1;
    
    zxca=0;

    end
    
    
    if rand>ruv
        
        num=num+1;
        
        veloC(num)=8;
     veloC(num)=veloC(num)*0.03;
    
    distance(num)=0;
 
    iter(num)=0;
    dragging(num) = 0;
    dr(num)=0;
    i(num)=1;
    aw(num)=0;
    xyz(num)=0;
    landing(num)=0;
    sv(num)=0;
    zv(num)=0;
    
    
    %% BALL POS
    
    if rand<0.25
        ballPos(num,:) = [round(200*rand) 4];
        angel(num)=rand*pi;
        xCord(num)=cos(angel(num))*veloC(num);
        yCord(num)=sin(angel(num))*veloC(num);
        
    elseif rand<0.5
        ballPos(num,:) = [round(200*rand) 96];
        angel(num)=-rand*pi;
        xCord(num)=cos(angel(num))*veloC(num);
        yCord(num)=sin(angel(num))*veloC(num);
        
    elseif rand<0.75
    ballPos(num,:) = [4 round(100*rand)];
    angel(num)=-rand*pi+pi/2;
        xCord(num)=cos(angel(num))*veloC(num);
        yCord(num)=sin(angel(num))*veloC(num);
        
    else
        ballPos(num,:) = [196 round(100*rand)];
        angel(num)=-rand*pi+pi/2;
        xCord(num)=cos(angel(num))*veloC(num);
        yCord(num)=sin(angel(num))*veloC(num);
        
    end
    [sv(num) zv(num)]=RoundOffAngle(angel(num));
    
    %%
    eval(['pathPos_' num2str(num) '(1,1:2)=ballPos(num,:)']);
    
    %% set PATH
    
    eval(['plot_' num2str(num) '=plot(pathPos_' num2str(num) '(:,1),pathPos_' num2str(num) '(:,2), ''lineStyle'', ''--'', ''linewidth'', 2)']);
                    
    
    %% Aeroplane
  
      
    h(num)=imagesc([ballPos(num,1)-ballRad ballPos(num,1)+ballRad], [ballPos(num,2)-ballRad ballPos(num,2)+ballRad], flipdim(sprite(sv(num)+1:sv(num)+100,zv(num)+1:zv(num)+100,:),3), 'alphadata', alphasprite(sv(num)+1:sv(num)+100,zv(num)+1:zv(num)+100,:));
    
    end
    
    
    
    
    
    ruv=ruv-0.0000001
    toc(start)
   pause(0.03-toc(start));
    
end