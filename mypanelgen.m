function [x_AF, z_AF] = mypanelgen(t_N,AOA,N)

t_N = input("The 4 digit NACA code is:" ,'s');
N = input("Number of Panels (Even Numbers Only):");
AOA = input('Angle of Attack is (Degrees):');
U_inf = input('Freestream Velocity is:');

%panelgen function that take NACA aerofoil code and the angle of attack to
%produce a plot of the chosen aerofoil
%Extract features from the NACAcode

% pulls the first digit out 
M_1 = str2double(t_N(1)); 
% pulls the second digit out 
P_2 = str2double(t_N(2));
% pulls the third and fourth digit out 
T_34 = str2double(t_N(3 : 4)); 

%Percentage values of airfoil properties
M = M_1/100;
P = P_2/10;
T = T_34/100;
 
%Number of grid points
if mod(N,2) == 1
    N = N + 1;
end
 
grid_points = N + 1;
 
%how finely to discretise distance from leading edge to trailing edge
for i = 1 : ceil(grid_points/2)
    x(i) = 1 - 0.5*(1 - cos(2*pi*((i - 1)/N)));
end
len = length(x);
 
%Camber function
y_c = zeros(1,len);
for i = 1 : len
    if (x(i) < P)
        y_c(i) = (M/P^2)*(2*P*x(i) - (x(i))^2);
 
    else
        y_c(i) = (M/((1-P)^2))*((1-2*P)+2*P*x(i) - (x(i))^2);
 
    end
end
 
%Thickness function
a_0 = 0.2969;
a_1 = -0.1260;
a_2 = -0.3516;
a_3 = 0.2843;
a_4 = -0.1036;
%a5 = -0.1015 doesn't allow for trailing and leading edge to connect

yt = zeros(1,len);
for j=1:len   
    yt(j) = 5*T*(a_0*sqrt(x(j)) + a_1*x(j) + a_2*(x(j))^2 + a_3*(x(j))^3 + a_4*((x(j)^4)));
end
 
%Coordinates
dycdx = zeros(1,len);
the_ta = zeros(1,len);
for r = 1 : len
    if x(r) < P
        dycdx(r) = (2*M/(P^2))*(P-x(r));
 
    else
        dycdx(r) = (2*M/((1 -(P^2))))*(P - (x(r)));
   
    end
    the_ta(r) = atan(dycdx(r));
       
end
x_u = zeros(1, len);
x_l = zeros(1, len);
z_u = zeros(1, len);
z_l = zeros(1, len);
for a=1:len
    x_u(a) = x(a) - yt(a)*sin((the_ta(a)));
    x_l(a) = x(a) + yt(a)*sin((the_ta(a)));
    z_u(a) = y_c(a) + yt(a)*cos((the_ta(a)));
    z_l(a) = y_c(a) - yt(a)*cos((the_ta(a)));
end
x_l = x_l(1 : end-1);
z_l = z_l(1 : end-1);
 
%Incorporating the Angle of Attack
magnitude_ = 99999999;
x_end = 99999999;
z_end = magnitude_ * tand(AOA) * -1;

x_AF = [x_l(1 : floor(grid_points/2)) flip(x_u(1 : ceil(grid_points/2))) 1 x_end];
z_AF = [z_l(1 : floor(grid_points/2)) flip(z_u(1 : ceil(grid_points/2))) 0 z_end];

 
%Plotting the airfoil
hold on
grid on
axis([-0.2 1.2 -0.7 0.7])
plot(x_AF(1 : end-2),z_AF(1 : end-2),'r-','LineWidth',3);
hold off
 