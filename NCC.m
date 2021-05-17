% ############################################################################################
% # Simulation of collective migration of neural crest cells. 
% # Author: Hamid Khataee.
% # Date: 01-December-2020.
% # Email: h.khataee@uq.edu.au
% # PLEASE DO NOT REDISTRIBUTE WITHOUT PERMISSION.
% #############################################################################################
clear all;
clc;
close all;

q = 10; % Flux rate.
R = 0.9; % interaction radius
eta = 1.5; % noise
v = 1; % cell velocity
Lx = 200; % domain length
Ly = 30;  % domain width
tsteps = 100000; % simulation steps.

dt = 1;
nparticles = 0;
c = [];                 % 1D binned density.
Cells_reached_end_of_channel = 0;
Start_reocrding_density = 0;
t_cells_filled_channel = 0;

theta = [];
xpos = [];
ypos = [];
va = [];
cross_time = [];
tstart = [];
ncross = 0;
tr = 1;
trace_1 = [];
ncontacts = [];

nparticles = nparticles+1; % number of particles
theta = [theta,2*pi*rand(1,1)]; % movement direction
xpos = [xpos,0*ones(1,1)]; % position_x of the first cell.
ypos = [ypos,Ly*rand(1,1)]; % position_y of the first cell.
tstart = [tstart,1*ones(1,1)]; % time that the first cell entered into the channel.

xpos_time = [];
ypos_time = [];
theta_time = [];
xpos_time(1,size(xpos,2)) = xpos;
ypos_time(1,size(ypos,2)) = ypos;
theta_time(1,size(theta,2)) = theta;

for i = 1:tsteps
    if floor (q) == q
        nparticles = nparticles+q; 
        theta = [theta,2*pi*rand(1,q)]; 
        xpos = [xpos,0*ones(1,q)]; 
        ypos = [ypos,Ly*rand(1,q)];
        tstart = [tstart,i*ones(1,q)]; 
    else
        nparticles = nparticles+floor(q);
        theta = [theta,2*pi*rand(1,floor(q))];
        xpos = [xpos,0*ones(1,floor(q))];
        ypos = [ypos,Ly*rand(1,floor(q))];
        tstart = [tstart,i*ones(1,floor(q))];
        if ~mod(i,round(1/(q - floor(q)))) % round q.
            nparticles = nparticles+1;
            theta = [theta,2*pi*rand(1,1)];
            xpos = [xpos,0*ones(1,1)];
            ypos = [ypos,Ly*rand(1,1)];
            tstart = [tstart,i*ones(1,1)];
        end
    end
        
    pos = [xpos.',ypos.'];
    idx = rangesearch(pos,pos,R);
    for j = randperm(nparticles)
        if xpos(j)<Lx % If cell j is withen the channel
            neighb = idx{j};            
            neighbtheta(j) = atan2(mean(sin(theta(neighb))),mean(cos(theta(neighb)))); 
            dtheta = eta*rand(1,1) - eta/2; 
             
            theta(j) = neighbtheta(j) + dtheta; 
            
            xpos(j) = xpos(j) + v*cos(theta(j))*dt; 
            ypos(j) = ypos(j) + v*sin(theta(j))*dt;
        end
    end
    
    for j = randperm(nparticles)     %Find list of particles that are about to hit boundary.
        if xpos(j) <0 
            theta(j) = pi-theta(j);
            xpos(j) = 0;
        end
        if ypos(j) <0
            ypos(j) = 0.1;
            theta(j) = 2*pi-theta(j);
        end
        if ypos(j) >Ly
            ypos(j) = Ly-0.1;
            theta(j) = 2*pi-theta(j);
        end
    end
    
    delVec = [];     %Delete Particles
    for j = randperm(nparticles)
        if xpos(j) > Lx
            ncross = ncross+1;
            cross_time(ncross) = i - tstart(j);
            delVec(end+1) = j;
            Cells_reached_end_of_channel = 1;           
        end
    end
    xpos(delVec) = [];
    ypos(delVec) = [];
    theta(delVec) = [];
    pos(delVec) = [];
    tstart(delVec) = [];
    
    if Cells_reached_end_of_channel == 1 && Start_reocrding_density == 0
        t_cells_filled_channel = i;        
        Start_reocrding_density = 1;
    end    
    xvel = [];     % Update velocity vectors.
    yvel = [];
    vel = [];
    xvel = v*cos(theta)*dt;
    yvel = v*sin(theta)*dt;
    vel = [xvel yvel];
    
    delVec_time = [];
    xvel_time = [];
    yvel_time = [];
    xpos_time(i,1:size(xpos,2)) = xpos;
    ypos_time(i,1:size(ypos,2)) = ypos;
    theta_time(i,1:size(theta,2)) = theta;
    delVec_time(i,1:size(delVec,2)) = delVec;
    xvel_time(i,1:size(xvel,2))= xvel;
    yvel_time(i,1:size(yvel,2))= yvel;
    
    bin_interval = 5;
    binedges = 0:bin_interval:Lx;
    c(i,:)= histcounts(xpos,binedges); 
    
    nparticles = size(theta,2);
    va(i) = (1/(nparticles*v))*norm(sum(vel));
    
    % Plots and Animations
    clf();
    subplot(2,1,1);
    hold on;
    if tr == 1 && nparticles>=1
        trace_1(i,:) = [xpos(1),ypos(1)];
        if xpos(1)> 0.99*Lx
            tr =0;
        end
    end
    quiver(xpos,ypos,2*xvel,2*yvel,0,'b');
    plot(trace_1(:,1),trace_1(:,2),'r','LineWidth',2);
    hold off;
    xlabel('x');
    ylabel('y');
    ax = gca;
    properties(ax);
    ax.TickLength = [0.03, 0.03]; 
    box on;     
    title(['v=',num2str(v),', R= ',num2str(R),', Eta= ',num2str(eta),', Flux= ',num2str(q),' (t = ',num2str(i),')']);
    axis([0,Lx,0,Ly]);
    subplot(2,1,2);
    plot(0.5:bin_interval:Lx,c(i,:),'ob');    
    legend(['t_{channel filled} = ', sprintf('%0.2f', t_cells_filled_channel)],'Location','northwest');
    xlabel('x');
    ylabel('Density c');
    ax = gca;
    properties(ax);
    ax.TickLength = [0.03, 0.03]; 
    box on;   
    pause(1/24);

figure('Renderer', 'painters', 'Position', [400 60 1125 150]);
hold on
xlabel('x');
ylabel('y');
quiver(xpos,ypos,2*xvel,2*yvel,0,'r');
ax = gca;
properties(ax);
axis([0,Lx,0,Ly]);
box on;
set(gca,'XTick',[], 'YTick', []);
hold off
close;
end

time = [1:tsteps]; 
figure;
hold on
xlabel('Time (MCS)'); 
ylabel('Order Parameter Va');
title(['v=',num2str(v),', R= ',num2str(R),', Eta= ',num2str(eta),', Flux= ',num2str(q),' (t = ',num2str(i),')']);
plot(time, va, 'ok','MarkerSize',12, 'linewidth', 1);
legend(['V_a = ', sprintf('%0.2f \\pm %0.2f', mean(va), std(va))],['All cells'], 'Location','southeast');
ax = gca;
properties(ax);
ax.TickLength = [0.03, 0.03]; 
ylim([0 1]);
set(gca,'YTick',[0 : 0.2 : 1]);
box on
hold off

figure;
hold on
xlabel('Number of cells exited from the channel'); 
ylabel('Time moved in the channel, Cross_{time}');
title(['v=',num2str(v),', R= ',num2str(R),', Eta= ',num2str(eta),', Flux= ',num2str(q),' (t = ',num2str(i),')']);
N_cells_exited_channel = [1:size(cross_time, 2)];
plot(N_cells_exited_channel, cross_time, 'or','MarkerSize',12, 'linewidth', 1);
legend(['Cross_{time} = ', sprintf('%0.2f \\pm %0.2f', mean(cross_time), std(cross_time))],['All cells'], 'Location','southeast');
ax = gca;
properties(ax);
ax.TickLength = [0.03, 0.03]; 
Nearest_divisible_byFive = round( max(cross_time) / 5 ) * 5;
box on
hold off

disp('Done!');
