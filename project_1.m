%% Project 1
clc;
clear all;
close all;
%% Mission Values
R=7370*6076; % nmi to ft
E=45/60; % min to hr
Alt = 45000; %ft
rho= 0.0003778; %slugs/ft^3
% Wing
Span = 212 + (7/12); %ft
Area = 4702; %ft^2
AR = Span.^2/Area;
e=0.7; % Osweld efficency
Cd_o = 0.023; % Parasitic Drag Coefficent
Ct = 0.55/3600; % 1/s
% Weights
W_Crew = 10*210; % lbs
W_Pass = 365*185; %lbs
W_Cargo = 365*40; %lbs
W_Payload = W_Crew + W_Pass +W_Cargo;

%% CL/CD and Velocity from Inital Weight Guesses
W_guess = linspace(300000,1000000,10); 
U_range = linspace(200,2000,20);
Cl_Cd_all = zeros(length(W_guess), length(U_range));  
Cl_Cd_max_vals = zeros(size(W_guess));          
U_max = zeros(size(W_guess));
for i = 1:length(W_guess)
    Wi = W_guess(i);
    Cl = Wi./(0.5*rho*(U_range.^2)*Area);
    Cd_i = (Cl.^2)./(AR*pi*e);
    Cl_Cd = Cl./(Cd_o + Cd_i);
    Cl_Cd_all(i,:) = Cl_Cd;
    [Cl_Cd_max_valls(i), idx] = max(Cl_Cd);
    U_max(i)= U_range(idx);
end


%% Plot CL/CD curves at each inital weight guess
figure;
hold on;
for i = 1:length(W_guess)
    plot(U_range,Cl_Cd_all(i,:),'Linewidth',1);
end
xlabel('Velocity (ft/s)');
ylabel('CL/CD');
title('CL/CD vs Velocity for Different Takeoff Weight Guesses');
legend(arrayfun(@(w) sprintf('W_{TO Guess} = %.0f lb', w), W_guess, 'UniformOutput', false), ...
       'Location', 'eastoutside', 'FontSize', 10);
grid on;


%% Weight Fractions and Iterations
W_converged = zeros(size(W_guess));
W1_W0 = 0.97;
W2_W1 = 0.985;
W5_W4 = 0.995;
for i = 1:length(W_guess)
    W_old = W_guess(i);
    iter = 0;
    while true
        iter = iter+1;
        [Cl_Cd_max, U_max] = find_Cl_Cd(W_old, rho, Area, AR, e, Cd_o);
        W3_W2 = exp((-1*R*Ct)./(Cl_Cd_max*U_max));
        W4_W3 = exp((-1*E*Ct)./(Cl_Cd_max));
        Wf_Wi = W1_W0*W2_W1*W3_W2*W4_W3*W5_W4;
        Wfuel_WMTOW = 1.06*(1-Wf_Wi);
        Wempty_WMTOW = 1.02*(W_old).^-0.06;
        WMTOW= W_Payload./(1-Wfuel_WMTOW-Wempty_WMTOW);
        if abs(WMTOW - W_old) < 1e3 || iter >= 50
            break;
        end
    W_converged(i)= WMTOW;
    end 
end
%% Plot Convergence
figure;
plot(W_guess, W_converged, 'b-o', 'LineWidth', 1.5); hold on;
plot(W_guess, W_guess, 'r--', 'LineWidth', 1.5);
xlabel('Guessed Takeoff Weight (lb)');
ylabel('Calculated Takeoff Weight (lb)');
title('Takeoff Weight Convergence');
legend('Calculated W_{TO}','Line of Equality','Location','best');
grid on;

%% Thrust
N_engines = 2; % 2 engines
W_conv = 767938; % lbs (converged value from iteration 5 based off convergence graph)
V = linspace(200,2000,400); % Velocity Range 
Drag_Force_parasitic = 0.5*rho*V.^2*Cd_o;
Drag_Force_induced = (2*W_conv.^2)./(rho*V.^2*Area*pi*e*AR);
Drag_Force_total = Drag_Force_induced + Drag_Force_parasitic;
Mach_numbers = [0.3 0.5 0.7 0.8 0.9 1.0];
Thrust_JT9D_1eng = [3500 3000 2700 2500 2300 2100];
Thrust_GE90_2eng = Thrust_JT9D_1eng * 1.92 * N_engines;
Temp = 216; % K (for air at 45000 ft)
R = 1716; % gas constant of air.
Gamma = 1.4; % gamma for air.
a = sqrt(Gamma*R*Temp); % speed of sound at 45000 ft.
V_Mach = a*Mach_numbers; % mach number to velocity conversion.
Thrust_available = interp1(V_Mach, Thrust_GE90_2eng, V, 'linear', 'extrap'); %finds velocity points in between mach conversion values.
[~, idx_int] = min(abs(Drag_Force_total - Thrust_available));
V_int = V(idx_int);
T_int = Drag_Force_total(idx_int);

figure; 
hold on; 
grid on;
plot(V, Drag_Force_parasitic, 'b--', 'LineWidth', 1.5);
plot(V, Drag_Force_induced, 'r--', 'LineWidth', 3);
plot(V, Drag_Force_total, 'k-', 'LineWidth', 1.5);
plot(V, Thrust_available, 'g-', 'LineWidth', 1.5);
plot(V_int, T_int, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);

xlabel('Velocity (ft/s)');
ylabel('Force (lbf)');
title('Thrust Required and Thrust Available at 45,000 ft');
legend('Parasitic Drag','Induced Drag','Total Drag (Thrust Required)', ...
       'Thrust Available (2 GE90 Engines)','Intersection Point','Location','best');
fprintf('Intersection point:\n');
fprintf('Velocity = %.1f ft/s (Mach %.2f)\n', V_int, V_int/a);
fprintf('Thrust Required = Thrust Available = %.0f lbf\n', T_int);

%% CL/CD function. 
function [Cl_Cd_max, U_max] = find_Cl_Cd(Wi, rho, Area, AR, e, Cd_o)
    U_range = linspace(200,2000,20);
    Cl = Wi./(0.5*rho*(U_range.^2)*Area);
    Cd_i = (Cl.^2)./(AR*pi*e);
    Cl_Cd = Cl./(Cd_o + Cd_i);
    [Cl_Cd_max, idx] = max(Cl_Cd);
    U_max= U_range(idx);
    end
