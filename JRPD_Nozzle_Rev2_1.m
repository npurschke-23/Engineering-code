clc
clear
close all

%inputs  
% ALL CALCULATIONS DONE IN SI UNITS
%this means that outputs are in meters
k = 2.22; 
F = 5000; %Desired Thrust (Newtons)
R = 311.64; %
T1 = 4000; %Temperature of combustion chamber (kelvin)
p1 = 2e6;  %chamber pressure (PA)
Rc = .0564; %radius of combustion chamber (meters)


p3 = 101325; %ambient pressure 
p2 = p3; %exit pressure -- for optimium expansion you make p2 = p3

%computing function defined in the bottom of the script
[AtA2,At,Ln,V2,Tt,mdot] = compute_nozzle(p1,p2,k,T1,R,F);


if numel(p1) > 1   %defined so that it only runs if we are looking at a variety of chamber pressures
%figures for data analysis of the nozzle at different chamber pressure
figure(1)
plot(p1,AtA2.^-1)
title("Area Ratio at Different Chamber Pressures")
ylabel("Area Ratio")
xlabel("Combustion Chamber Pressure (PA)")

figure(2)
plot(p1,10000*At)
title("Throat Area at Different Chamber Pressures")
ylabel("Throat Area (cm^2)")
xlabel("Combustion Chamber Pressure (PA)")

figure(3)
plot(p1,100*Ln)
title("Length of Nozzle at Different Chamber Pressures")
ylabel("Length of Nozzle (cm)")
xlabel("Combustion Chamber Pressure (PA)")

figure(4)
plot(p1,V2)
title("Exit Velocity at Different Chamber Pressures")
ylabel("Exit Velocity (m/s)")
xlabel("Combustion Chamber Pressure (PA)")

figure(5)
plot(p1,mdot)
title("Mass Flow Rates at Different Chamber Pressures")
ylabel("Mass Flow Rate (kg/s)")
xlabel("Combustion Chamber Pressure (PA)")

figure(6)
plot(p1,10000.*At./AtA2)
title("Exit Area at Different Chamber Pressures")
ylabel("Exit Area (cm^2)")
xlabel("Combustion Chamber Pressure (PA)")

%figure(7)
%plot(p1,rt)
%title("Throat Radii at Different Chamber Pressures")
%label("Throat Radius (cm)")
%xlabel("Combustion Chamber Pressure (PA)")
end



%this section is utlized to get a cad file as well as a graph of the
%resulting geometry of the nozzle

Rt = sqrt(At./pi());    %throat radius
Lcn = Rc+(1.5*(sqrt(2)-1)-1).*Rt;  %length from chamber to throat
Ae = At/AtA2;     %exit area
Re = sqrt(Ae./pi()); %exit radius

%this loop only runs when a single chamber pressure is selected
if isscalar(p1) 


% use figure 3-13 to solve for throat angle towards exit and exit angle
%based on area ratio 5 i = 32.5  e = 13.5
%these are required for the selected bell curve
theta3 = deg2rad(32.5); % theta in radians
theta4 = deg2rad(13.5); % theta in radians

%defining the bounds as we change from different equations to solve for the
%radius (y) from the centerline of the nozzle along its length (x)
x0 = -(Rc+(1.5*(sqrt(2)-1)-1)*Rt);
x1 = -(1.5*sqrt(2)/2)*Rt;
x2 = 0;
x3 = .45*Rt*cos((pi()/2)-theta3);
y3 = 1.45*Rt-sqrt((.45*Rt)^2-(x3)^2); 
x4 = Ln;
y4 = Re; 

%from chamber to starting curve of throat
xs01 = linspace(x0, x1, 300);
ys01 = -(xs01-x0)+Rc;

%from starting curve of throat to throat
xs12 = linspace(x1, x2, 300);
ys12 = 2.5*Rt-sqrt((1.5*Rt).^2-xs12.^2);

%from throat to ending curve of throat
xs23 = linspace(x2, x3, 300);
ys23 = 1.45*Rt-sqrt((.45*Rt).^2-xs23.^2);

%from ending curve of throat to exit
%this utilized a 3rd degreen polynomial to define the geometry of the bell
%curve of the nozzle selected

%notably the two theta values ensure tangency to the selected angles at the
%entrance and exit of this section
m1 = tan(theta3);
m2 = tan(theta4);

% A*c = b where c = [a3; a2; a1; a0] Matrice Setup handcalculated
A = [ x3^3, x3^2, x3, 1;
      x4^3, x4^2, x4, 1;
      3*x3^2, 2*x3, 1, 0;
      3*x4^2, 2*x4, 1, 0 ];

b = [ y3; y4; m1; m2 ];

c = A\b;   % solves for coefficients [a3; a2; a1; a0]

%defining the function in terms of x and y for plotting
xs34 = linspace(x3, x4, 300);
ys34 = polyval(c, xs34);


%plotting the 4 different segments of the nozzle into 1 
figure(6)
hold on  

plot(xs01, ys01, 'r', 'LineWidth', 1.5) 
plot(xs12, ys12, 'b', 'LineWidth', 1.5)  
plot(xs23, ys23, 'g', 'LineWidth', 1.5) 
plot(xs34, ys34, 'k', 'LineWidth', 1.5)  
yline(0, '--k', 'LineWidth', 1);   %centerline for better visualization

hold off   %adjust the bounds if large changes for better plot
xlim([-.05 .07])
ylim([-.02 .1])
xlabel('length (m)')
ylabel('radius (m)')
title('Nozzle Geometry')

%combining all the nozzle points to make a asc file with 3d points. This
%will be saved in the same file as this live script and can be imported
%into cad, fit curve and revolve for a model of the nose cone
x = [xs01,xs12,xs23,xs34];
y = [ys01,ys12,ys23,ys34];
z = zeros(1200, 1);
Nozzle_Geometry = [x.',y.',z];

save('Nozzle_Geometry.asc', "Nozzle_Geometry", '-ascii')

end


%outputs this is all from page 60-61 in Sutton Textbook
% computes:
%   nozzle area ratio
%   nozzle throat area
%   length of nozzle
%   velocity at exit
%   temperature at throat
%   mass flow rate
function [AtA2,At,Ln,V2,Tt,mdot] = compute_nozzle(p1,p2,k,T1,R,F)
p_ratio = p2./p1;   %pressure ratio
AtA2 = (((k+1)/2).^(1/(k-1))).*(p_ratio.^(1/k)).*sqrt(((k+1)./(k-1)).*(1-(p_ratio).^((k-1)./k))); %inverse of nozzle area ratio
Tt = (2*T1)/(k+1);  %throat temperature
V2 = sqrt((2*k*R*T1/(k-1)).*(1-p_ratio.^((k-1)/k))); %exit velocity
mdot = F./V2;   %mass flow rate, same at all cross sections
At = (mdot./p1).*sqrt((R*T1)/(k*(2/(k+1)).^((k+1)/(k-1)))); %throat area
rt = sqrt(At/pi());            %from area to radius  %throat radius
%Ln = .8*(((sqrt(AtA2.^-1))-1).*rt)./tand(15);  %bell curve length equation for 80% length bell curve
re = sqrt((rt.^2)./AtA2);
Ln = (0.8).*((re-rt)./tand(15)); %bell curve length equation for 80% length bell curve
end