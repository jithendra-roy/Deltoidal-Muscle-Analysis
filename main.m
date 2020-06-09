clc;
clear all;
height = input("Enter height of individual(in m): ");
weight = input("Enter weight of individual(in kg): ");
max_force_angle = input("Enter maximum arm angle(in deg): ");
external_weight = input("Enter external load(in kg): ");
Lf = height * 0.156;    % length of forearm
Lu = height * 0.174;    % length of upper arm
Mf = weight * 0.016;    % weight of forearm
Mu = weight * 0.028;    % weight of upper arm
arm_length = Lf + Lu;   % length of total arm
arm_weight = Mu + Mu;     % weight of total arm
[Fjx_list, Fjy_list, T0_list, Pcom_list, sigma_list] = calc(max_force_angle, arm_length, external_weight, arm_weight, Lf, Lu, Mf, Mu);
x_data = 1:max_force_angle;
% COM radial distance
figure
plot(x_data,Pcom_list);
title('Arm angle(Deg) vs COM radial distance')
xlabel('Arm Angle')
ylabel('Pcom')
% Deltoidal Muscle Force
figure
plot(x_data,T0_list);
title('Arm angle(Deg) vs Deltoidal muscle tension(N)')
xlabel('Arm Angle')
ylabel('T0')
% Deltoidal Muscle Stress
figure
plot(x_data, sigma_list);
title('Arm angle(Deg) vs Muscle stress(N/m^2)')
xlabel('Arm Angle')
ylabel('sigma')
% Shoulder Joint Stress Y
figure
plot(x_data,Fjy_list);
title('Arm angle(Deg) vs Shoulder Y joint Force(N)')
xlabel('Arm Angle')
ylabel('Fjy')
% shoulder Joint Stress X
figure
plot(x_data,Fjx_list);
title('Arm angle(Deg) vs Shoulder X joint force(N)')
xlabel('Arm Angle')
ylabel('Fjx')

function [Fjx_list, Fjy_list, T0_list, Pcom_list, sigma_list] = calc(max_force_angle, arm_length, external_weight, arm_weight, Lf, Lu, Mf, Mu)
    Lm = 0.0645;
    g = 9.81; % acceleration due to gravity
    t = 1:(max_force_angle);
    nt = length(t);
    Fjx_list = zeros(nt,1);
    Fjy_list = zeros(nt,1);
    T0_list = zeros(nt,1);
    % T1_list = zeros(nt,1);
    sigma_list = zeros(nt,1);
    Pcom_list = zeros(nt,1);
    delt_angle = linspace(0,10,max_force_angle);
    for angle=1:max_force_angle
        % Define where positions are
        arm_pos = (Lf * Mf + Lu * Mu) / arm_weight;
        external_pos = arm_length;
        % COM weight
        %Wcom = external_weight + arm_weight;
        % COM pos changes with arm rotation
        Pcom = (arm_weight * arm_pos + external_weight * external_pos)/(arm_weight+external_weight) * sin(angle);
        PCSA0 = 0.0017;     % For Deltoids
        syms T0 sig Fjy Fjx
        % Stress Eqn
        eqn_stress_deltoids = T0-PCSA0*sig == 0;
        % Force in Y
        eqn_y = cos(angle-delt_angle(angle))*T0-arm_weight*g-external_weight*g-Fjy == 0;
        % Force in X
        eqn_x = Fjx - sin(angle-delt_angle(angle))*T0 == 0;
        % Moment Equation
        eqn_m = sin(delt_angle(angle))*Lm*T0-arm_weight*g*arm_pos*sin(angle)-g*external_weight*external_pos*sin(angle) == 0;
        [A,B] = equationsToMatrix([eqn_y, eqn_x, eqn_m, eqn_stress_deltoids], [Fjx, Fjy, T0, sig]);
        output = linsolve(A,B);
        temp = output;
        Fjx_list(angle) = temp(1,1);
        Fjy_list(angle) = temp(2,1);
        T0_list(angle) = temp(3,1);
        sigma_list(angle) = temp(4,1);
        Pcom_list(angle) = Pcom;
    end
end