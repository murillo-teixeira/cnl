clear

%% Parameters 

p1 = -0.021151;
p2 = 0.092551;
p3 = -0.014188;
p4 = 0.077947;

%% Equilibrium

d_eq = 0.005;
u_eq = 0;

x1_eq = -p2*d_eq/p1;
x2_eq = -(u_eq + p4)/(p3 - x1_eq);

%% Phase portrait

odefun = @(t, x) [p1*x(1) + p2*d_eq;(p3 - x(1))*x(2) + u_eq + p4];
plotpp(odefun)

%% Running the model
% 
% d = 0;
% u = 0;
% 
% simout = sim("simulation.slx");
% tout = simout.tout;
% logsout = simout.logsout;
% 
% x1 = logsout.get('x1').Values.Data;
% x2 = logsout.get('x2').Values.Data;
% 
% plot(tout, x1, tout, x2)