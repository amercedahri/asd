APPENDIX. MATLAB CODE
VEHICLE_DYNAMICS
% REMUS parameters developed by Chris Churan edited by Tyler Furukawa
clear
clc
z_g = 1.96e-2; x_b = 0; W = 299; buoy = 306;
I_z = 3.45; I_y = 3.45; I_x = 1.77e-1;
U = 1.5; to = 0; tf = 80;
m = 299/9.81; M_q = -6.87; M_qdot = -4.88;
M_w = 30.7; M_wdot = -1.93; M_d = -34.6;
Z_q = -9.67; Z_qdot = -1.93;
Z_w = -66.6; Z_wdot = -35.5; Z_d = -50.6;
% Dynamics ------------------------------------------------------------
% modified for [wr; q; theta]
M = [m-Z_wdot -Z_qdot 0; -M_wdot I_y-M_qdot 0; 0 0 1];
A_0 = [Z_w m*U+Z_q 0; M_w M_q -z_g*W; 0 1 0];
B_0 = [Z_d; M_d; 0];
A = inv(M)*A_0
B = inv(M)*B_0
C = [0 0 1]
%Check open loop poles
open_loop_poles = eig(A)
%Check Controlability
Control=[B,A*B,A^2*B];
Controllable=rank(Control)
%Check Observability
Observe=[C',A'*C',A'^2*C'];
Observability=rank(Observe)
p = [-.5+.866i, -.5-.866i, -1] %using butterworth pattern
K = place(A,B,p) %[-1.6807, 0.2935, 0.0889]
TRACKING
function command=tracking(x,z)
aim = 4; %aims 4 meters ahead
%%%%%Normal Altitude Command%%%%%%%%%
Z_mine=35;
min_clear=2;
X_mine = 60;
sigma=5; %sigma^2=25
47
%Gaussian Function
x_g1 = x; %x-position
z_g1 = Z_mine-min_clear*exp(-(X_mine-x_g1)^2/(2*sigma^2));% z-posit
x_g2 = x+0.1; %creates a delta x of 0.1 meter
z_g2 = Z_mine-min_clear*exp(-(X_mine - x_g2)^2/(2*sigma^2));
%Tangent
T2=[(x_g2-x_g1);(z_g2-z_g1)];
%Normal
N2=[-(z_g2-z_g1);(x_g2-x_g1)];
%Position
P2=[(0);(z-z_g1)];
%cross track error
CTE=P2'*N2/sqrt(N2'*N2);
command=-atan2(-CTE,aim);
TRACKING_WITH_SLOPE
function command=tracking_with_slope(x,z)
aim = 4; %aims 4 meters ahead
Z_mine=35;
min_clear=2;
X_mine = 60;
sigma=5; %sigma^2=25
%Gaussian Function
x_g1 = x; %x-position
z_g1 = Z_mine-min_clear*exp(-(X_mine-x_g1)^2/(2*sigma^2)); % z-posit
x_g2 = x+0.1; %creates a delta x of 0.1 meter for slope calculation
z_g2 = Z_mine - min_clear*exp(-(X_mine - x_g2)^2/(2*sigma^2));
%Tangent
T2=[(x_g2-x_g1);(z_g2-z_g1)];
%Normal
N2=[-(z_g2-z_g1);(x_g2-x_g1)];
%Position
P2=[(0);(z-z_g1)];
%cross track error
CTE=P2'*N2/sqrt(N2'*N2);
slope=atan2(-(z_g2-z_g1),(x_g2-x_g1));
command=slope-atan2(-CTE,aim);
48
TRACKING_LOOKAHEAD
function command=tracking_lookahead(x,z)
aim = 4; %aims 4 meters ahead
Z_mine=35;
min_clear=2;
X_mine = 60;
sigma=5; %sigma^2=25
lookahead = 4.58 %ideal lookahead = 4.58 m
%Gaussian Function
x_g1 = x+lookahead; %modified x-position
z_g1 = Z_mine-min_clear*exp(-(X_mine-x_g1)^2/(2*sigma^2)); %mod z-posit
x_g2 = x_g1+0.1; %creates a delta x of 0.1 meter
z_g2 = Z_mine - min_clear*exp(-(X_mine - x_g2)^2/(2*sigma^2));
%Tangent
T2=[(x_g2-x_g1);(z_g2-z_g1)];
%Normal
N2=[-(z_g2-z_g1);(x_g2-x_g1)];
%Position
P2=[(0);(z-z_g1)];
%cross track error
CTE=P2'*N2/sqrt(N2'*N2);
command=-atan2(-CTE,aim);
CONTROLLER ERRORS
%calculates and plots controller errors
x_sim=(1:1:120);
[num_points, columns] = size(x_pos);
E=zeros(num_points,1);
for k=1:1:num_points
E(k,1)=sqrt(((ocean_depth-z_pos(k))-(Mine_altitude+min_clear*exp(-
(x_pos(k)-X_mine)^2/(2*sigma^2))))^2);
end
figure(1);clf
plot(x_pos,E), title('Vertical Error')
xlabel('Horizontal Position (m)')
ylabel('Vertical Error (m)')
%Calculates and plots sum of error
sum_error=zeros(num_points,1);
for l=2:1:num_points
sum_error(l,1)=sum_error(l-1)+E(l-1);
end
figure(2);
49
plot(x_pos,sum_error), title('Total Vertical Error')
xlabel('Horizontal Position (m)')
ylabel('Total Vertical Error (m)')
max_error = max(E)
avg_error = sum(E)/num_points
total_error =sum(E)
OCEAN MODELS
Wall_Z = 5; %height of obstacle
Wall_X = 100;
Wall_Z2 = 4.5;
Wall_X2= 110;
%Floor 1
X_Floor1 = 0:Interval:Wall_X-Interval;
Z_Floor1 = zeros(1,length(X_Floor1));
%Floor 2
Z_Floor2 = 0:Interval:Wall_Z;
X_Floor2 = Wall_X*ones(1,length(Z_Floor2));
%Floor 3
X_Floor3 = Wall_X+Interval:Interval:2*Wall_X;
Z_Floor3 = zeros(1,length(X_Floor3));
%Floor 4
Z_Floor4 = 0:Interval:Wall_Z2;
X_Floor4 =Wall_X2*ones(1,length(Z_Floor4));
%plots the floor
subplot(2,1,1)
plot(X_Floor1,Z_Floor1,'g*'), hold on
plot(X_Floor2,Z_Floor2,'g*')
plot(X_Floor3,Z_Floor3,'g*')
plot(X_Floor4,Z_Floor4,'g*')
axis ([0 120 0 60])
SONAR_MODEL
Wall_Z = 5; %height of obstacle
Wall_X = 60;
Interval = 0.1;
%Floor 1
X_Floor1 = 0:Interval:Wall_X-Interval;
Z_Floor1 = zeros(1,length(X_Floor1));
%Floor 2
Z_Floor2 = 0:Interval:Wall_Z;
X_Floor2 = Wall_X*ones(1,length(Z_Floor2));
%Floor 3
X_Floor3 = Wall_X+Interval:Interval:120;
Z_Floor3 = zeros(1,length(X_Floor3));
50
%plots the floor
subplot(2,1,1)
plot(X_Floor1,Z_Floor1,'g*'), hold on
plot(X_Floor2,Z_Floor2,'g*')
plot(X_Floor3,Z_Floor3,'g*')
axis ([0 120 0 60])
title('Sonar Model')
%plots ordered path
x_sim=(1:1:120);
for j=1:1:length(x_sim)
plot(x_sim(j),Mine_altitude+min_clear*exp(-(x_sim(j)
X_mine)^2/(2*sigma^2)),'r*');hold on
end
X_Remus = 0; Z_Remus = 5; %Remus Position (Can take from SIMULINK)
plot(X_Remus, Z_Remus,'bo')
Theta = 0; %Remus Pitch in degrees(taken from SIMULINK model)
Max_Sonar_Range = 100; %in meters
%Calculates Slopes of Sonar Beam
sonar_angle =-22.5+Theta:-.5:-45+Theta; %interval determines # of beams
for a=1:length(sonar_angle)
sonar_slope(a)=tand(sonar_angle(a));
%plots sonar lines
Beam_X=[X_Remus, X_Remus+Max_Sonar_Range*cosd(-sonar_angle(a))];
Beam_Z=[Z_Remus, Z_Remus-Max_Sonar_Range*sind(-sonar_angle(a))];
plot(Beam_X,Beam_Z,'k:')
end
%finds the highest beam that intercepts the wall
for b=1:length(sonar_angle)
Wall_Intercept(b)=sonar_slope(b)*(Wall_X-X_Remus) + Z_Remus;
if Wall_Intercept(b)<=Wall_Z & Wall_Intercept(b)>=0
Beam_Intercept = b;
Highest_Beam_Z = Wall_Intercept(b);
break
els Highest_Beam_Z=0;e
end
end
%find occlusion area
if Highest_Beam_Z>0
Occlusion_X=[Wall_X, Wall_X, Wall_X, Wall_X+(Max_Sonar_Range-
(Wall_X-X_Remus))*cosd(-sonar_angle(Beam_Intercept))];
Occlusion_Z=[0, Highest_Beam_Z, Highest_Beam_Z, Highest_Beam_Z-
(Max_Sonar_Range-(Wall_X-X_Remus))*sind(-sonar_angle(Beam_Intercept))];
else Occlusion_X = [0, 0]; Occlusion_Z=[0, 0];
end
subplot(2,1,2)
plot(X_Floor1,Z_Floor1,'g*'), hold on
plot(X_Floor2,Z_Floor2,'g*')
plot(X_Floor3,Z_Floor3,'g*')
plot(Occlusion_X,Occlusion_Z,'r--')
axis ([0 120 0 60])
51
title('Occlusion Area')
%y=mx+b for the sonar b = Z_remus and Global_X = x +X_Remus
%Find where bottom beam intercepts the bottom 0=mx+Z_remus or x = -
Z_remus/m
Lowest_Beam_X = -Z_Remus/sonar_slope(length(sonar_angle));
if sqrt(Z_Remus^2 +Lowest_Beam_X^2) < Max_Sonar_Range &
Lowest_Beam_X+X_Remus > X_Remus %check range & if x intercept < X_Remus
Left_X = Lowest_Beam_X +X_Remus;
els Left_X = 0;e
end
%if bottom beam doesn't intercept bottom, find its “Z” at wall
if Left_X<=0 & abs(Wall_X-X_Remus)<=Max_Sonar_Range
%y = mx+b where x = (Wall_X-X_Remus), and b = Z_Remus
Lowest_Beam_Z=sonar_slope(length(sonar_angle))*(Wall_X-X_Remus) +
Z_Remus;
if Lowest_Beam_Z<=Wall_Z;
Low_Wall_Intercept = Lowest_Beam_Z;
else
Low_Wall_Intercept = 0;
end
end
%if X intercept is behind wall, find its “Z” at wall
if Left_X>Wall_X & abs(Wall_X-X_Remus)<=Max_Sonar_Range
%y = mx+b where x = Wall_X, and b = Z_Remus
Lowest_Beam_Z=sonar_slope(length(sonar_angle))*(Wall_X-X_Remus) +
Z_Remus;
if Lowest_Beam_Z<=Wall_Z;
Low_Wall_Intercept = Lowest_Beam_Z;
else
Low_Wall_Intercept = 0;
end
end
SONAR_MOVING
clf
Wall_X = 100; Wall_Z = 5; %height of obstacle
Wall_X2= 110; Wall_Z2 = 4.5;
Interval = 0.1;
Max_Sonar_Range = 100; %in meters
Gaus2_offset = 35;
Gaus2_range = Wall_X-Gaus2_offset;
Orig_alt = 3;
%Floor 1
X_Floor1 = 0:Interval:Wall_X-Interval;
Z_Floor1 = zeros(1,length(X_Floor1));
%Floor 2
Z_Floor2 = 0:Interval:Wall_Z;
52
X_Floor2 = Wall_X*ones(1,length(Z_Floor2));
%Floor 3
X_Floor3 = Wall_X+Interval:Interval:2*Wall_X;
Z_Floor3 = zeros(1,length(X_Floor3));
%Floor 4
Z_Floor4 = 0:Interval:Wall_Z2;
X_Floor4 =Wall_X2*ones(1,length(Z_Floor4));
x_sim=(0:.2:2*Wall_X);
figure(1),clf;
for k=1:1:length(x_sim)
if x_sim(k)<Gaus2_range
%plot orig alt
plot(x_sim(k),Orig_alt,'r*');hold on
else if x_sim(k)<(Wall_X+Gaus2_offset)
%plot gaussian
plot(x_sim(k),Orig_alt+(min_clear+Wall_Z-Orig_alt)*exp(-
(x_sim(k)-Wall_X)^2/(2*sigma2^2)),'r*')
else
plot(x_sim(k),Orig_alt,'r*')
end
end
end
%plot(X_mine,Mine_altitude,'ro');
plot(x_pos,(ocean_depth-z_pos));
axis([0 2*Wall_X 0 15]), title('Ordered versus Actual')
xlabel('X - meters'), ylabel('Altitude - meters')
mov = avifile('sonar_Movie.avi','Compression','Cinepak','FPS',1);
for c=3:2:length(x_pos);
set(gcf,'doublebuffer','on');
figure(c)
%plots the floor
subplot(2,1,1)
plot(X_Floor1,Z_Floor1,'g*'), hold on
plot(X_Floor2,Z_Floor2,'g*')
plot(X_Floor3,Z_Floor3,'g*')
plot(X_Floor4,Z_Floor4,'g*')
%plots actual course
%plot(x_pos,ocean_depth-z_pos,'c*')
axis ([30 130 0 15])
title('Sonar Model')
X_Remus = x_pos(c);
Z_Remus = ocean_depth-z_pos(c); %Remus Posit (Taken from SIMULINK)
plot(X_Remus, Z_Remus,'bo')
Theta = Remus_theta(c); %Remus Pitch in degrees(taken from SIMULINK
%Calculates Slopes of Sonar Beam
53
sonar_angle =(-2.5+Theta):-Interval:(-25+Theta); %determines # of
beams
for a=1:length(sonar_angle)
sonar_slope(a)=tand(sonar_angle(a));
%plots sonar lines
Beam_X=[X_Remus, X_Remus+Max_Sonar_Range*cosd(-
sonar_angle(a))];
Beam_Z=[Z_Remus, Z_Remus+Max_Sonar_Range*sind(sonar_angle(a))];
plot(Beam_X,Beam_Z,'k:')
end
%finds the highest beam that intercepts the wall
for b=1:length(sonar_angle)
Wall_Intercept(b)=sonar_slope(b)*(Wall_X-X_Remus) + Z_Remus;
if Wall_Intercept(b)<=Wall_Z & Wall_Intercept(b)>=0
Beam_Intercept = b;
Highest_Beam_Z = Wall_Intercept(b);
break
else
Beam_Intercept=-1;
Highest_Beam_Z=-1;
end
end
%find if “sees” second obstacle
%if intercepts 1st wall
Sonar_Wall_Z2=-1;%initialize each time
Sonar_Wall_X2=-1;%intiialize each time
if Beam_Intercept>1;
Wall2_Intercept = sonar_slope(Beam_Intercept)*(Wall_X2-X_Remus)
+Z_Remus;
if Wall2_Intercept <= Wall_Z2
Sonar_Wall_Z2=Wall2_Intercept:Interval:Wall_Z2;
Sonar_Wall_X2=Wall_X2*ones(1,length(Sonar_Wall_Z2));
end
elseif Beam_Intercept<1;
%check to see if intercept 2nd wall
for d=length(sonar_angle):-1:1 %cycle from lowest beam
Wall_Intercept2=sonar_slope(d)*(Wall_X2-X_Remus) + Z_Remus;
if Wall_Intercept2<=Wall_Z2 & Wall_Intercept2>=0 &
X_Remus<Wall_X2
%Beam_Intercept2 = d;
%Lowest_Beam_Z2 = Wall_Intercept2(d);
Sonar_Wall_Z2=Wall_Intercept2:Interval:Wall_Z2;
Sonar_Wall_X2=Wall_X2*ones(1,length(Sonar_Wall_Z2));
eakbr
else
Sonar_Wall_X2=-1;
Sonar_Wall_Z2=-1;
end
end
end
if Highest_Beam_Z>0 & X_Remus<Wall_X
%find occlusion area
54
Occlusion_X=[Wall_X, Wall_X+(Max_Sonar_Range-(Wall_X-
X_Remus))*cosd(-sonar_angle(Beam_Intercept))];
Occlusion_Z=[Highest_Beam_Z, Highest_Beam_Z-(Max_Sonar_Range-
(Wall_X-X_Remus))*sind(-sonar_angle(Beam_Intercept))];
els Occlusion_X = [0, 0]; Occlusion_Z=[0, 0];e
end
subplot(2,1,2)
plot(Occlusion_X,Occlusion_Z,'r--'), hold on
axis ([30 130 0 15])
title('Sonar Image w/Occlusion Area')
%y=mx+b for the sonar b = Z_remus and Global_X = x +X_Remus
%Find where top beam intercepts the bottom 0=mx+Z_remus or x = -
Z_remus/m
Top_Beam_X=-Z_Remus/sonar_slope(1);
if sqrt(Z_Remus^2+Top_Beam_X^2)<Max_Sonar_Range &
Top_Beam_X+X_Remus<Wall_X & Top_Beam_X+X_Remus>X_Remus%check range and
if x intercept < Wall_X
Right_X=Top_Beam_X+X_Remus;
else Right_X=-1;
end
%Find where bottom beam intercepts the bottom 0=mx+Z_remus or x = -
Z_remus/m
Lowest_Beam_X = -Z_Remus/sonar_slope(length(sonar_angle));
if sqrt(Z_Remus^2+Lowest_Beam_X^2)<Max_Sonar_Range &
Lowest_Beam_X+X_Remus<Wall_X & Lowest_Beam_X+X_Remus>X_Remus%check
range and if x intercept<Wall_X
Left_X = Lowest_Beam_X +X_Remus;
if Right_X>0
X_Sonar_floor1=Left_X:Interval:Right_X;
else
X_Sonar_floor1=Left_X:Interval:Wall_X-Interval;
end
Z_Sonar_floor1=zeros(1,length(X_Sonar_floor1));
Z_Sonar_floor2=0:Interval:Highest_Beam_Z;
X_Sonar_floor2=Wall_X*ones(1,length(Z_Sonar_floor2));
elseif abs(Wall_X-X_Remus)<=Max_Sonar_Range
%if bottom beam doesn't intercept bottom or is behind wall,
find its “Z” at wall
Lowest_Beam_Z=sonar_slope(length(sonar_angle))*(Wall_X-X_Remus)
+ Z_Remus;
if Lowest_Beam_Z<=Wall_Z;
Low_Wall_Intercept = Lowest_Beam_Z;
Z_Sonar_floor1=Low_Wall_Intercept:Interval:Highest_Beam_Z;
X_Sonar_floor1=Wall_X*ones(1,length(Z_Sonar_floor1));
Z_Sonar_floor2=-1;
X_Sonar_floor2=-1;
else Z_Sonar_floor1 = -1; X_Sonar_floor1=-1; Z_Sonar_floor2=-1;
X_Sonar_floor2=-1;
end
end
%subplot(3,1,3)
plot(X_Sonar_floor1,Z_Sonar_floor1,'g*')
plot(X_Sonar_floor2,Z_Sonar_floor2,'g*')
55
plot(Sonar_Wall_X2,Sonar_Wall_Z2,'g*')
axis ([30 130 0 15])
F = getframe(gcf);
mov = addframe(mov,F);
pause(.1);
end
mov = close(mov);
TRACKING_POPUP
function command=tracking_popup(x,z)
aim = 4; %aims 4 meters ahead
X_mine = 100;
Z_mine=35; %Altitude of 5 (when ocean_depth=40)
org_depth = 37; %altitude of 3 (when ocean_depth=40)
%for popup Gaussian
ho=5;
sigma1=5;
%for Obstacle Avoidance Gaussian
min_clear=3;
sigma2=12;
clear = (min_clear+org_depth-Z_mine);
gaus1_offset=80;
gaus1_range=X_mine-gaus1_offset;
gaus2_offset = 40; %do gaussian 40m before and after obstacle
gaus2_range = X_mine-gaus2_offset;
lookahead = 4.5; %”look” ahead 4.5
delta_x = .1; %for slope calculation
if x<(gaus1_range-lookahead)
%do original altitude
x_g1 = x+lookahead; %modified x-position
z_g1 = org_depth;
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
z_g2 = org_depth;
elseif x<(gaus2_range-lookahead)
%do 1st gaus
x_g1 = x+lookahead; %modified x-position
z_g1 = org_depth - ho*exp(-(gaus1_range+25-x_g1)^2/(2*sigma1^2));
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
z_g2 = org_depth - ho*exp(-(gaus1_range+25-x_g2)^2/(2*sigma1^2));
elseif x<(X_mine+gaus2_offset)
%do gaussian
x_g1 = x+lookahead; %modified x-position
z_g1 = org_depth - clear*exp(-(X_mine-x_g1)^2/(2*sigma2^2));
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
56
z_g2 = org_depth - clear*exp(-(X_mine - x_g2)^2/(2*sigma2^2));
else
x_g1 = x+lookahead; %modified x-position
z_g1 = org_depth;
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
z_g2 = org_depth;
end
%Tangent
T2=[(x_g2-x_g1);(z_g2-z_g1)];
%Normal
N2=[-(z_g2-z_g1);(x_g2-x_g1)];
%Position
P2=[(0);(z-z_g1)];
%cross track error
CTE=P2'*N2/sqrt(N2'*N2);
command=-atan2(-CTE,aim);
TVD CALCULATION
Original_Depth=37; %altitude =3m m
Deviation=zeros(1,length(x_pos));
Area=Deviation;
for n=2:1:length(x_pos);
Deviation(n)=(Original_Depth - (z_pos(n)+z_pos(n-1))/2); %avg
deviation between 2 pts
Area(n)=Deviation(n)*(x_pos(n)-x_pos(n-1)); %tapezoid rule
end
Max_Deviation=max(Deviation)
TVD=sum(Area)
TRACKING SPLINE
function out=tracking_spline2(x,z)
aim = 4; %aims 4 meters ahead
X_mine = 100;
Z_mine=35; %Altitude of 5 (when ocean_depth=40)
org_depth = 37; %altitude of 3 (when ocean_depth=40)
%for modified popup
ho=.75;
sigma1=2;
%for obstacle avoidance Gaussian
min_clear=2;
sigma2=5;
clearance = org_depth-(Z_mine-min_clear);
gaus1_offset=80;
gaus1_range=X_mine-gaus1_offset;
57
gaus2_offset = 20; %do gaussian 20m before and after obstacle
gaus2_range = X_mine-gaus2_offset;
lookahead = 4.58; %”look” ahead 4.58
delta_x = .1; %for slope calculation
%second obstacle
Wall_Z2 = 35.5;
Wall_X2= 110;
clear2=org_depth-(Wall_Z2-min_clear);
Orig_alt = 3;
%for spline
Interval=0.1;
spline_offset=.1;
spline_range=X_mine+spline_offset;
X_slope1=spline_range+lookahead;
spline_X_int=Wall_X2-X_slope1;
X_spline1=X_slope1;
X_spline2=X_spline1+(spline_X_int)/3;
X_spline3=X_spline2+(spline_X_int)/3;
X_spline4=Wall_X2;
X_spline5=Wall_X2+10;
X_spline6=X_spline5+1;
X_spline7=X_spline6+1;
X_spline8=2*X_mine;
X_spline=[X_spline1 X_spline2 X_spline3 X_spline4 X_spline5 X_spline6
X_spline7 X_spline8];
Z_spline1=clearance*exp(-(X_mine-X_slope1)^2/(2*sigma2^2))*(100-
X_slope1)/-(sigma2^2); %slope of spline
Z_spline2=org_depth-clearance*exp(-(X_mine-
X_slope1)^2/(2*sigma2^2));%height at gaussian
Z_spline3=Z_spline2-(Z_spline2-(Wall_Z2-min_clear))/3;
Z_spline4=Z_spline3-(Z_spline2-(Wall_Z2-min_clear))/3;
Z_spline5=Wall_Z2-min_clear;
Z_spline6=org_depth;
Z_spline7=Z_spline6;
Z_spline8=Z_spline6;
Z_spline9=Z_spline6;
Z_spline10=0.0;
Z_spline=[Z_spline1 Z_spline2 Z_spline3 Z_spline4 Z_spline5 Z_spline6
Z_spline7 Z_spline8 Z_spline9 Z_spline10];
%xx=X_slope1:Interval:X_spline8;
%yy=spline(X_spline,Z_spline,xx);
if x<(gaus1_range-lookahead)
%do original altitude
x_g1 = x+lookahead; %modified x-position
z_g1 = org_depth;
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
z_g2 = org_depth;
elseif x<(gaus2_range-lookahead)
%do 1st gaus
58
x_g1 = x+lookahead; %modified x-position
z_g1 = org_depth - ho*exp(-(gaus1_range+25-x_g1)^2/(2*sigma1^2));
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
z_g2 = org_depth - ho*exp(-(gaus1_range+25-x_g2)^2/(2*sigma1^2));
elseif x<(spline_range)
%do 2nd gaussian
x_g1 = x+lookahead; %modified x-position
z_g1 = org_depth - clearance*exp(-(X_mine-x_g1)^2/(2*sigma2^2));
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
z_g2 = org_depth - clearance*exp(-(X_mine - x_g2)^2/(2*sigma2^2));
else
x_g1 = x+lookahead; %modified x-position
z_g1 = spline(X_spline,Z_spline,x_g1);
x_g2 = x_g1+delta_x; %creates a delta x of 0.1 meter
z_g2 = spline(X_spline,Z_spline,x_g2);
end
%Tangent
T2=[(x_g2-x_g1);(z_g2-z_g1)];
%Normal
N2=[-(z_g2-z_g1);(x_g2-x_g1)];
%Position
P2=[(0);(z-z_g1)];
%cross track error
error=P2'*N2/sqrt(N2'*N2);
slope=atan2((z_g2-z_g1),(x_g2-x_g1));
command= -atan2(-error,aim);
x_g1=x_g1;
z_g1=z_g1;
out=[command,x_g1,z_g1];
