%% Inhibitor Model Simulations and Graph Plotting Code

clear all;
clc;

%%% Parameters

% Host Cell Area, and Exudate Height and Volume
p.h = 10^(-1);             % cm    Height
p.A_r = 49;                % cm^2  Area
p.V = p.h*p.A_r;           % cm^3  Volume

% Number of Binding Sites Occupied by a Bacterium or an Inhibitor
p.phi_Bac = 1;             % sites cell^(-1)
p.phi_A   = 1;             % sites inhibitor^(-1)

% Hill Coefficient
p.n = 1;                   % dimensionless

%%% Initial Conditions

B_Free_init  = 5*10^6/p.V;       % bacteria cm^(-3)    Free Bacteria
B_Bound_init = 0;                % bacteria cm^(-2)    Bound Bacteria
A_Free_init  = 3*10^8/p.V;       % inhibitors cm^(-3)  Free Inhibitors
A_Bound_init = 0;                % inhibitors cm^(-2)  Bound Inhibitors
p.E_init     = 1.261*10^8/p.A_r; % sites cm^(-2)       Binding Sites

% Initial Condition Vectors for ODE Solvers
IC_vec         = [B_Free_init,B_Bound_init];
IC_vec_Treat   = [B_Free_init,B_Bound_init,A_Free_init,A_Bound_init];

%%% Load and Set Fitted Parameters

load Fitted_Parameters;

Set = 2;

% Case A: Sets 1-2  - Set 2 is presented in the main text.
% Case B: Sets 3-7  - Set 3 is presented in the main text.
% Case C: Sets 8-11 - Set 8 is presented in the main text.
% Case D: Set  12   - Set 12 is presented in the main text.

p.r_F       = Fitted_Parameter_Table(1,Set);  % hr^(-1)
p.r_B       = Fitted_Parameter_Table(2,Set);  % hr^(-1)
p.K_F       = Fitted_Parameter_Table(3,Set);  % cells cm^(-3)
p.K_B       = Fitted_Parameter_Table(4,Set);  % cells cm^(-2)
p.alpha_Bac = Fitted_Parameter_Table(5,Set);  % hr^(-1) sites^(-1)
p.beta_Bac  = Fitted_Parameter_Table(6,Set);  % hr^(-1)
p.delta_B   = Fitted_Parameter_Table(7,Set);  % hr^(-1)
p.eta_max   = Fitted_Parameter_Table(8,Set);  % dimensionless
p.gamma     = Fitted_Parameter_Table(9,Set);  % sites cm^(-2)
p.psi_Bac   = Fitted_Parameter_Table(10,Set); % hr^(-1)
p.alpha_A   = Fitted_Parameter_Table(11,Set); % hr^(-1) sites^(-1)
p.beta_A    = Fitted_Parameter_Table(12,Set); % hr^(-1)
p.psi_A     = Fitted_Parameter_Table(13,Set); % hr^(-1)

%%% Set time parameters and vectors over which to solve the ODEs

% Long Time Span
t_end_lomg       = 3000;              % hr
t_int_long       = 0.1;               % hr
t_grid_full_long = 0:t_int_long:3000; % hr

%%% Simulations

options = odeset('AbsTol',1e-20,'RelTol',1e-10,'MaxOrder',2); % To avoid negative soln (free bac)

% No Treatment

[~,sol] = ode15s(@(t,y)No_Treatment_ODE_Fn(t,y,p),t_grid_full_long,IC_vec);

B_Free_sol_NT_long   = sol(:,1)*p.V;
B_Bound_sol_NT_long  = sol(:,2)*p.A_r;
E_sol_NT_long        = p.E_init*p.A_r - p.phi_Bac*B_Bound_sol_NT_long;
B_T_NT_long          = B_Free_sol_NT_long + B_Bound_sol_NT_long;

% Single Inhibitor Dose

[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn(t,y,p),t_grid_full_long,IC_vec_Treat);
        
B_Free_sol_Treat_long   = sol(:,1)*p.V;
B_Bound_sol_Treat_long  = sol(:,2)*p.A_r;
A_Free_sol_Treat_long   = sol(:,3)*p.V;
A_Bound_sol_Treat_long  = sol(:,4)*p.A_r;
E_sol_Treat_long        = p.E_init*p.A_r - p.phi_Bac*B_Bound_sol_Treat_long - p.phi_A*A_Bound_sol_Treat_long;
B_T_Treat_long          = B_Free_sol_Treat_long + B_Bound_sol_Treat_long;


% Debridement Only

% New time parameters and vectors are required as the solver must be
% interrupted upon debridement
t_end_D     = 24;                         % hr
t_int_D     = 0.1;                        % hr
t_grid_48   = 0:t_int_D:2*t_end_D;        % hr
t_grid_loop = 0:t_int_D:t_end_D;          % hr
length_t_grid_loop = length(t_grid_loop); % hr
days = 100;                               % hr
t_grid_full_D = 0:t_int_D:t_end_D*days;   % hr

IC_vec_loop = IC_vec;

sol_full = zeros(length_t_grid_loop*days-(days-1),2);

[~,sol] = ode15s(@(t,y)No_Treatment_ODE_Fn(t,y,p),t_grid_48,IC_vec_loop,options);

sol_full(1:length_t_grid_loop*2 -  1,:) = sol;

IC_vec_loop = [0,sol(end,2)];

for i = 3:days
    
[~,sol] = ode15s(@(t,y)No_Treatment_ODE_Fn(t,y,p),t_grid_loop,IC_vec_loop,options);

sol_full(length_t_grid_loop*(i-1)+(2-i):length_t_grid_loop*i-(i-1),:) = sol;

IC_vec_loop = [0,sol(end,2)];

end

B_Free_D_NT   = sol_full(:,1)*p.V;
B_Bound_D_NT  = sol_full(:,2)*p.A_r;
E_num_D_NT    = p.E_init*p.A_r - p.phi_Bac*B_Bound_D_NT;

B_Total_D_NT = B_Free_D_NT + B_Bound_D_NT;


% Debridement with a Single Inhibitor Dose 

IC_vec_loop_Treat = IC_vec_Treat;

sol_full = zeros(length_t_grid_loop*days-(days-1),4);

[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn(t,y,p),t_grid_48,IC_vec_loop_Treat,options);

sol_full(1:length_t_grid_loop*2 -  1,:) = sol;

IC_vec_loop_Treat = [0,sol(end,2),0,sol(end,4)];

for i = 3:days
    
[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn(t,y,p),t_grid_loop,IC_vec_loop_Treat,options);

sol_full(length_t_grid_loop*(i-1)+(2-i):length_t_grid_loop*i-(i-1),:) = sol;

IC_vec_loop_Treat = [0,sol(end,2),0,sol(end,4)];

end

B_Free_D_Treat   = sol_full(:,1)*p.V;
B_Bound_D_Treat  = sol_full(:,2)*p.A_r;
A_Free_D_Treat   = sol_full(:,3)*p.V;
A_Bound_D_Treat  = sol_full(:,4)*p.A_r;
E_num_D_Treat    = p.E_init*p.A_r - p.phi_Bac*B_Bound_D_Treat - p.phi_A*A_Bound_D_Treat;

B_Total_D_Treat = B_Free_D_Treat + B_Bound_D_Treat;

% Debridement with Regular Inhibitor Doses

IC_vec_loop_Treat = IC_vec_Treat;

sol_full = zeros(length_t_grid_loop*days-(days-1),4);

[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn(t,y,p),t_grid_48,IC_vec_loop_Treat,options);

sol_full(1:length_t_grid_loop*2 -  1,:) = sol;

IC_vec_loop_Treat = [0,sol(end,2),IC_vec_Treat(3),sol(end,4)];

for i = 3:days
    
[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn(t,y,p),t_grid_loop,IC_vec_loop_Treat,options);

sol_full(length_t_grid_loop*(i-1)+(2-i):length_t_grid_loop*i-(i-1),:) = sol;

IC_vec_loop_Treat = [0,sol(end,2),IC_vec_Treat(3),sol(end,4)];

end

B_Free_D_Treat_Reg   = sol_full(:,1)*p.V;
B_Bound_D_Treat_Reg  = sol_full(:,2)*p.A_r;
A_Free_D_Treat_Reg   = sol_full(:,3)*p.V;
A_Bound_D_Treat_Reg  = sol_full(:,4)*p.A_r;
E_num_D_Treat_Reg    = p.E_init*p.A_r - p.phi_Bac*B_Bound_D_Treat_Reg - p.phi_A*A_Bound_D_Treat_Reg;

B_Total_D_Treat_Reg = B_Free_D_Treat_Reg + B_Bound_D_Treat_Reg;

% Regular Inhibitor Doses Only

IC_vec_loop_Treat = IC_vec_Treat;

sol_full = zeros(length_t_grid_loop*days-(days-1),4);

[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn(t,y,p),t_grid_48,IC_vec_loop_Treat,options);

sol_full(1:length_t_grid_loop*2 -  1,:) = sol;

IC_vec_loop_Treat = [sol(end,1),sol(end,2),sol(end,3)+IC_vec_Treat(3),sol(end,4)];

for i = 3:days
    
[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn_NC(t,y,p),t_grid_loop,IC_vec_loop_Treat,options);

sol_full(length_t_grid_loop*(i-1)+(2-i):length_t_grid_loop*i-(i-1),:) = sol;

IC_vec_loop_Treat = [sol(end,1),sol(end,2),sol(end,3)+IC_vec_Treat(3),sol(end,4)];

end

B_Free_ND_Treat_Reg_NC   = sol_full(:,1)*p.V;
B_Bound_ND_Treat_Reg_NC  = sol_full(:,2)*p.A_r;
A_Free_ND_Treat_Reg_NC   = sol_full(:,3)*p.V;
A_Bound_ND_Treat_Reg_NC  = sol_full(:,4)*p.A_r;
E_num_ND_Treat_Reg_NC    = p.E_init*p.A_r - p.phi_Bac*B_Bound_ND_Treat_Reg_NC - p.phi_A*A_Bound_ND_Treat_Reg_NC;

B_Total_ND_Treat_Reg_NC = B_Free_ND_Treat_Reg_NC + B_Bound_ND_Treat_Reg_NC;

%%% Constant Debridement, with and without MAM7 Treatment

% High Clearance Parameters for Constant Debridement
p.psi_Bac_2 = 1000;
p.psi_A_2 = 1000;

% Constant Debridement & No Inhibitor Treatment

sol_full = zeros(length(t_grid_full_D),2);

t_grid_remaining = 0:t_int_D:t_end_D*(days-2);

[~,sol] = ode15s(@(t,y)No_Treatment_ODE_Fn(t,y,p),t_grid_48,IC_vec);


sol_full(1:length(t_grid_48),:) = sol;

IC_vec_2 = [0,sol(end,2)];

[~,sol] = ode15s(@(t,y)No_Treatment_Const_Deb_ODE_Fn(t,y,p),t_grid_remaining,IC_vec_2);

sol_full(length(t_grid_48):end,:) = sol;

B_Free_sol_CD_NT   = sol_full(:,1)*p.V;
B_Bound_sol_CD_NT  = sol_full(:,2)*p.A_r;
E_sol_CD_NT        = p.E_init*p.A_r - p.phi_Bac*B_Bound_sol_CD_NT;
B_T_CD_NT          = B_Free_sol_CD_NT + B_Bound_sol_CD_NT;

% Constant Debridement & Single Inhibitor Dose

sol_full = zeros(length(t_grid_full_D),4);

t_grid_remaining = 0:t_int_D:t_end_D*(days-2);


[~,sol] = ode15s(@(t,y)Inhib_Dose_ODE_Fn(t,y,p),t_grid_48,IC_vec_Treat);

sol_full(1:length(t_grid_48),:) = sol;

IC_vec_2 = [0,sol(end,2),0,sol(end,4)];

[T,sol] = ode15s(@(t,y)Inhib_Dose_Const_Deb_ODE_Fn(t,y,p),t_grid_remaining,IC_vec_2);
        
sol_full(length(t_grid_48):end,:) = sol;

B_Free_sol_CD_Treat   = sol_full(:,1)*p.V;
B_Bound_sol_CD_Treat  = sol_full(:,2)*p.A_r;
A_Free_sol_CD_Treat   = sol_full(:,3)*p.V;
A_Bound_sol_CD_Treat  = sol_full(:,4)*p.A_r;
E_sol_CD_Treat        = p.E_init*p.A_r - p.phi_Bac*B_Bound_sol_CD_Treat - p.phi_A*A_Bound_sol_CD_Treat;
B_T_CD_Treat          = B_Free_sol_CD_Treat + B_Bound_sol_CD_Treat;

%% Graph Plotting

%%% Comparing NT and Treat

figure;
plot(t_grid_full_long,B_T_NT_long,'k--'); hold on;
plot(t_grid_full_long,B_T_Treat_long,'r');
xlabel('time, t (hr)');
ylabel(['total num. bacteria,',sprintf('\n'),' B_T (cells)']);
legend('No Treatment','Single Inhib. Dose');
axis([0 3000 0 1.5*10^8]);

%%% Comparing Dependent Variables

% No Treatment

% Single Inhibitor Dose
figure;

subplot(1,2,1);
plot(t_grid_full_long,B_Free_sol_NT_long,'k--'); hold on;
plot(t_grid_full_long,B_Bound_sol_NT_long,'b');
plot(t_grid_full_long,E_sol_NT_long,'m:');
title('No Treatment');
xlabel('time, t (hr)');
ylabel('num. cells/sites');
xlim([0 1500]);
ylim([0 1.5*10^8]);

subplot(1,2,2);
plot(t_grid_full_long,B_Free_sol_Treat_long,'k--'); hold on;
plot(t_grid_full_long,B_Bound_sol_Treat_long,'b');
plot(t_grid_full_long,A_Free_sol_Treat_long,'r-.');
plot(t_grid_full_long,A_Bound_sol_Treat_long,'g');
plot(t_grid_full_long,E_sol_Treat_long,'m:');
title('Single Inhib. Dose');
xlabel('time, t (hr)');
ylabel(['num.',sprintf('\n'),'cells/inhib./sites']);
legend('B_F','B_B','A_F','A_B','E');
xlim([0 1500]);
ylim([0 1.5*10^8]);

%%% Comparing Treatment Strategies

figure;

subplot(1,2,1);
plot(t_grid_full_long,B_T_NT_long,'k--'); hold on; % t_grid_full_D B_T_NT_medium
plot(t_grid_full_D,B_Total_D_NT,'b-.');
plot(t_grid_full_long,B_T_Treat_long,'r'); % t_grid_full_D B_T_Treat_medium
plot(t_grid_full_D,B_Total_D_Treat,'g');
plot(t_grid_full_D,B_Total_ND_Treat_Reg_NC,'m:');
plot(t_grid_full_D,B_Total_D_Treat_Reg,'c');
xlabel('time, t (hr)');
ylabel(['total num. bacteria,',sprintf('\n'),' B_T (cells)']);
legend('No Treatment','Reg. Deb.','Inhib.',...
    'Inhib. & Reg. Deb.','Reg. Inhib.','Reg. Inhib. & Reg. Deb.');
xlim([0 2000]);
ylim([0 1.5*10^8]);

subplot(1,2,2);
plot(t_grid_full_long,B_T_NT_long,'k--'); hold on; % t_grid_full_D B_T_NT_medium
plot(t_grid_full_D,B_T_CD_NT,'b-.');
plot(t_grid_full_long,B_T_Treat_long,'r'); % t_grid_full_D B_T_Treat_medium
plot(t_grid_full_D,B_T_CD_Treat,'g');
xlabel('time, t (hr)');
ylabel(['total num. bacteria,',sprintf('\n'),' B_T (cells)']);
legend('No Treatment','Cont. Deb.','Inhib.','Inhib. & Cont. Deb.');
xlim([0 2000]);
ylim([0 1.5*10^8]);

