%% Microscopic Model

clear;close all;clc;

myVideo = VideoWriter('robot'); %open video file
% myVideo.FrameRate = 30;  %can adjust this
open(myVideo)

%%

k       = 1;
lengthh = 4;
width   = 4;
sep     = 1.0;
r       = lengthh/sep;
c       = width/sep;
n       = r*c;

Pg{n}  = [];

grid_initial_x(1,1) = 0.5;
grid_initial_z(1,1) = 0.5;

for j=1:r
    for i = 1:c
        Pg{k}               = [grid_initial_x(i);grid_initial_z(j)];
        grid_initial_x(i+1) = grid_initial_x(i)+sep;
        grid_initial_z(j+1) = grid_initial_z(j)+sep;
        k = k+1;
    end
end

x = (cell2mat(Pg))';

nRobots = 5;
colors = {'r','b','k','g','m'};

figure(1);hold on; h = gscatter(x(:,1),x(:,2));

d = imread('House.png');
image(flipud(d), 'XData', [1.75 2.25], 'YData', [4 5]);

ylim([0;5.5]);

xx = 1.75-0.5/4;
yy = 4.0;

for i = 1:nRobots
    
    xx = xx+0.5/4;
    Robots{1,i} = [xx yy];
    
end

m = 2;

while ~(isempty(x))
    
    for i = 1:nRobots
        
        if isempty(x); break; end
        
        [~,dist] = dsearchn(Robots{m-1,i},x);
        weights  = 1./dist;
        idx      = randsample(1:size(x,1),1,'true',weights);
        Robots{m,i} = x(idx,:); x(idx,:)  = [];
        
    end
    
    m = m+1;
end

for i = 1:nRobots
    
    Robot{i}   = cell2mat(Robots(:,i));
    
%     Works in MATLAB 2020 onli
%     xdata      = interp(Robot{i}(:,1),1.5);
%     ydata      = interp(Robot{i}(:,2),1.5);
%     Robot{i}   = [xdata ydata];
    
    datalen{i} = 1:size(Robot{i},1);
    
    % plot(Robot{i}(:,1),Robot{i}(:,2),'-o','MarkerSize',8);
    
end

for i = 1:nRobots
   for j = 1:size(Robot{i},1)-1
       x1{j,1} = linspace(Robot{i}(j,1),Robot{i}(j+1,1),20)';
       y1{j,1} = linspace(Robot{i}(j,2),Robot{i}(j+1,2),20)';
   end
   x2 = cell2mat(x1);
   y2 = cell2mat(y1);
   Robot2{i} = [x2 y2];
   datalen{i} = 1:size(Robot2{i},1);
end

Robot = Robot2;
pause(5.0);
for k = 1:numel(datalen{1})
    
    for i = 1:nRobots
        try
            p(i) = plot(Robot{i}(1:k,1),Robot{i}(1:k,2),...
                'color',colors{i},'LineStyle','-',...
                'Marker','*','MarkerSize',6,...
                'LineWidth',2);
        catch
            continue;
        end
    end
    pause(0.3);
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
end

legend([p(1) p(2) p(3) p(4) p(5)],{'Robot 1','Robot 2','Robot 3','Robot 4','Robot 5'});

xlabel('x');
ylabel('y');
title('Robots'' position at t = 3');

close(myVideo)