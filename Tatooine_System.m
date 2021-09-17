%% Tatooine System
%
%  written by Mikey Lambert
%  December 5th, 2015
%
%  Analytical Mechanics Project

%%
clear

G = 6.667E-11;      % gravitational constant
MS1 = 1.7751E30;    % mass of Tatoo I
MS2 = 2.2051E30;    % mass of Tatoo II
for n=1:1000
    RS1(n) = 6.507E8/1000;  % radius of each sun in km
    RS2(n) = 6.995E8/1000;
end
% convert to kilometers
c1 = 1.8122E9/1000;      % c values for solar orbits
c2 = 1.4588E9/1000;
us = 3.1323E29;          % reduced solar mass
% solar eccentricities
es1 = 0.5;
% formula on pg. 5
es2 = sqrt((-c2/c1)*(1-es1^2)+1);

AU = 150E6; % Astromical unit in km
theta = linspace(0,2*pi,1000);  % polar independent variable

%% Polar Solar Graph

% orbital figures for Tatooine
aT = 1.62E11/1000;              % semi-major axis in km
eT = 0.5;                       % eccentricity (arbitrary)
cT = aT*(1-eT^2);               % c-value for Tatooine
rT = cT./(1+eT.*cos(theta));    % distance from barycenter

% calculate solar orbital periods (both are equal)
% use meters, not kilometers for these calculations
TS = sqrt(4*pi^2*(c1*1000/(1-es1^2))^3/(G*(MS1+MS2)));
% compare to period of Tatooine
TT = sqrt((4*pi^2*(aT*1000)^3)/(G*(MS1+MS2)));
z = TT/TS;

% attempt to keep suns on a separate variable from Tatooine
Stheta = linspace(0,2*pi*z,1000);  

r1 = c1./(1+es1.*cos(theta));       % solar distances from barycenter
r2 = -c2./(1+es2.*cos(theta));

polar(theta,r1)
hold on
polar(theta,r2)
title('Tatoo Binary Star Orbits in km')
text(-1.25*c1,-1.25*c1,'Tatoo I')
text(c1,-1.2*c1,'Tatoo II')
text(0,-1.5*c1,['Solar Orbital Period: ' num2str(TS/3600) ' hrs.'])
hold off

%% Solar System Animation

% animate
figure
t = 0;
x1 = r1.*cos(Stheta);       % cartesians for Tatoo I
y1 = r1.*sin(Stheta);
x2 = r2.*cos(Stheta);       % cartesians for Tatoo II
y2 = r2.*sin(Stheta);
xT = rT.*cos(theta);        % cartesians for Tatooine
yT = rT.*sin(theta);

    for n = 1:1000
        plot(x1(n),y1(n),'o')
        hold on
        plot(x2(n),y2(n),'o')
        plot(xT(n),yT(n),'o')
        hold off
        title('Tatooine Orbital System')
        text(-0.35*AU,-0.85*AU,['Distance to CM of Binary Star: ' num2str((sqrt(xT(n)^2+yT(n)^2))/AU) ' AU'])
        text(-0.2*AU,-1.025*AU,['Tatooine Orbital Period: ' num2str(TT/86400) ' days'])
        text(-0.02*AU,-1.2*AU,'*planetry/solar radii are not to scale')
        xlim([-2.5E8 2.5E8])
        xlabel('km')
        ylim([-2E8 2E8])
        ylabel('km')
        
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [A,cm] = rgb2ind(im,256);
        if n == 1;
            imwrite(A,cm,'TatooineSystem.gif','gif', 'Loopcount',inf);
        else
            imwrite(A,cm,'TatooineSystem.gif','gif','WriteMode','append');
        end
    n = n+1;
    end


%% Barycenter Animation

figure
t = 0;
x1 = r1.*cos(theta);
x2 = r2.*cos(theta);
y1 = r1.*sin(theta);
y2 = r2.*sin(theta);
while t == 0
    for n = 1:100
        plot(x1(10*n),y1(10*n),'o')
        hold on
        plot(x2(10*n),y2(10*n),'o')
        hold off
        xlim([-5E6 5E6])
        ylim([-5E6 5E6])
        title('Tatoo I and II')
        text(0.55*c1,-2.5*c1,'*solar radii are not to scale')
        text(0.15*c1,-2*c1,['Solar Orbital Period: ' num2str(TS/3600) ' hrs.'])
        drawnow
    end
end

%%

while t == 0
    for n = 1:500
        P = polar(theta, 2.5E8*ones(size(theta)))
        set(P,'Visible','off')
        hold on
        polar(Stheta(n*2), r1(n*2), 'o')
        polar(Stheta(n*2), r2(n*2), 'o')
        polar(theta(n*2), rT(n*2), 'o')
        hold off
        title('Tatooine Orbital System')
        text(-0.35*AU,-0.75*AU,['Distance to CM of Binary Star: ' num2str((rT(n*2)/AU)) ' AU'])
        text(-0.2*AU,-1.125*AU,['Tatooine Orbital Period: ' num2str(TT/86400) ' days'])
        text(-0.02*AU,-1.3*AU,'*planetry/solar radii are not to scale')
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [A,cm] = rgb2ind(im,256);
        if n == 1;
            imwrite(A,cm,'TatooineSystem.gif','gif', 'Loopcount',inf);
        else
            imwrite(A,cm,'TatooineSystem.gif','gif','WriteMode','append');
        end
    n = n+1;
    end
end
