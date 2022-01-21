
%% Numerical simulation of 1-D compressible flow
%% Study Case: Sod Shock Tube
%% Type: Approximate Riemann solution
%% Solve the 1-D Euler equation
%% Using Compressible Solver 
%% Involve flux splitting (FVS / FDS) + high-resolution upwind schemes (UPW) / shock capturing schemes (SCS)
%% Version: 2.0
%% Date: 2021/12/21

%% HW: Fundamentals of Computational Fluid Dynamics
%% Name: Feng Zhenghao
%% College: College of Engineering
%% ID: 2101112008

%% Program Initiation
clear
clf
close all
clc

warning off

%% Global Variables Definition

% global L Gamma R Cv;

% Characteristic length of the tube (x = [-L/2, L/2])
L = 1.0;

% Flow Paremeters
Gamma = 1.4;
R = 286.9;  % J/kg * K
% R = 8.314;
Cv = R / (Gamma - 1);            % Specific heat at constant volume
Cp = (Gamma * R) / (Gamma - 1);  % Specific heat at constant pressure

%% Grid Generation
N = 201;  % Number of grids space in the x coordinate
xp = linspace(-L / 2, L / 2, N)';  % x coordinate of grid points
dx = L / (N - 1);

xp_mid = ceil(N / 2); % The grid number of x = 0 position

% Arrays of Properties
u_arr = zeros(N, 1);     % velocity
rho_arr = zeros(N, 1);   % density
p_arr = zeros(N, 1);     % pressure

rho = zeros(N, 1);
u = zeros(N, 1);
p = zeros(N, 1);
E = zeros(N, 1);
T = zeros(N, 1);
c = zeros(N, 1);

% Time interval (values can be adjusted to ensure the stability)
dt = 0.001;

% Max. cal. round
max_step = 100;

% Max. elapsed time
max_tot_time = max_step * dt;

%% Pre-processing: I/O File Setting
% Set the global I/O file folder path
savefolder = ['Program_Sod_Shock_Tube', '_MaxTime_',num2str(max_tot_time, '%.3f')];
save_output_folder = ['.\', savefolder, '\'];
mkdir(save_output_folder);

%% Setting Initial Condition
% when x <  0, (ul, rhol, pl) = (0.0, 1.0  , 1.0);
% when x >= 0, (ur, rhor, pr) = (0.0, 0.125, 0.1);
u_arr(1 : (xp_mid-1)) =   0.0;
rho_arr(1 : (xp_mid-1)) = 1.0;
p_arr(1 : (xp_mid-1)) =   1.0;

u_arr(xp_mid : end) =   0.0;
rho_arr(xp_mid : end) = 0.125;
p_arr(xp_mid : end) =   0.1;

%% Pre-processing: Selected Setting & Algorithm
% Specify the Flux Splitting Method
% 1 - Flux Vector Splitting (FVS)
% 2 - Flux Differnce Splitting (FDS)
flag_flu_spl = 1;

    % Family of FVS methods
    % 1 - Steger-Warming (S-W)
    % 2 - Lax-Friedrich (L-F)
    % 3 - Van Leer
    % 4 - Liou-Steffen Splitting - Advection Upstream Splitting (AUSM) Method
    flag_fvs_met = 1;
    
    % Specify the Flux Reconstruction Method
    % 1 - Direct Reconstruction (for F(U))
    % 2 - Characteristic Reconstruction
    flag_flu_rec = 1;
    
    % Family of FDS methods
    % FDS - Roe scheme
    flag_fds_met = 1;

% Specify the Algorithm of High-Resoluation Flux Spatial discretization (Shock Capturing Scheme)
% For Complete Flux Reconstruction
% Avilable types
% 1 - General Upwind & Compact Schemes (forward / backward)
% 2 - Special Shock-Capturing Schemes
flag_spa_typ = 1;

    % Global
    % [General Upwind & Compact Schemes (forward / backward)]
    % 1 - 1_od upwind scheme (2 points) [Safest! Treated as Standard Benchmark]
    % 2 - 2_od upwind scheme (3 points)
    % 3 - 3_od upwind scheme (4 points with bias)
    % 4 - 5_od upwind scheme (6 points with bias)
    % compact scheme (five points with 2_od C.S.)
    flag_upw_typ = 1;
    
    % [Special Shock-Capturing Schemes]
    % Godunov
    % 1 - TVD - Total Variation Diminishing Scheme (Van Leer Limiter)
    % 2 - NND - Non-oscillatory, Non-free-paremeters Dissipative Difference Scheme
    % 3 - WENO - Weighted Essentially Non-Oscillatory Method (e.g. 5 od WENO proposed by Jiang & Shu)
    % MUSCL - Monotone Upstream-Centered Schemes for Conservation Laws
    flag_scs_typ = 1;

% Specify the Algorithm for temporal marching [OK]
% Global
% 1 - Euler time marching
% 2 - Trapezoid formula
% 3 - 2_od Runge-Kunta method (Heun formula)
% 4 - 3_od TVD Runge-Kunta method
% 5 - 4_od Runge-Kunta method
flag_tim_mar = 4;

% Zone title set (Export results)
zone_title_FVM = ["FVS", "FDS"];
    zone_title_FVS = ["S-W", "L-F", "Van Leer", "AUSM"];
    zone_title_FDS = ["Roe"];
zone_title_SPA = ["UPW", "SCS"];
    zone_title_UPW = ["1_od", "2_od", "3_od", "5_od"];
    zone_title_SCS = ["TVD (VL Limiter)", "NND", "WENO (5_od)"];
zone_title_MAR = ["Euler", "Trape", "R-K (2_od)", "R-K (3_od TVD)", "R-K (4_od)"];

if (flag_flu_spl == 1)
    zone_title_FVM_comb = [char(zone_title_FVM(flag_flu_spl)), '_', char(zone_title_FVS(flag_fvs_met))];
elseif (flag_flu_spl == 2)
    zone_title_FVM_comb = [char(zone_title_FVM(flag_flu_spl)), '_', char(zone_title_FDS(flag_fds_met))];
end
if (flag_spa_typ == 1)
    zone_title_SPA_comb = [char(zone_title_SPA(flag_spa_typ)), '_', char(zone_title_UPW(flag_upw_typ))];
elseif (flag_spa_typ == 2)
    zone_title_SPA_comb = [char(zone_title_SPA(flag_spa_typ)), '_', char(zone_title_SCS(flag_scs_typ))];
end
zone_title_MAR_comb = char(zone_title_MAR(flag_tim_mar));

zone_title_comb = [zone_title_FVM_comb, '_', zone_title_SPA_comb, '_', zone_title_MAR_comb];

% Export movie settings
flag_exp_mov = 0;
flag_exp_avi = 0;
flag_exp_gif = 1;

%% Array Initiation
lambda = zeros(3, 1);     % eigenvalues matrix 
lambda_p = zeros(3, 1);  % pos. eigenvalues matrix
lambda_n = zeros(3, 1);  % neg. eigenvalues matrix

U = zeros(N, 1);     % variables vector
% F = zeros(N, 1);     % invisid vector
F_p = zeros(N, 3);   % pos. invisid vector
F_n = zeros(N, 3);   % neg. invisid vector
Fx = zeros(N, 3);    % invisid term
Fx_p = zeros(N, 3);  % pos. invisid term
Fx_n = zeros(N, 3);  % neg. invisid term

Fh_p = zeros((N - 1), 3); % half point (j + 1/2) invisid vector
Fh_n = zeros((N - 1), 3); % half point (j + 1/2) invisid vector

% Cal. the initial rho, u, p, T, c
for i = 1 : N

    % Step 1: Cal. rho, u, p, T, c
    rho = rho_arr(i);
    u = u_arr(i);
    p = p_arr(i);

    T = p / (rho * R);
    c = sqrt(Gamma * p / rho);

    E = rho * ((Cv * T) + (0.5 * u * u));

    U(i, 1) = rho;
    U(i, 2) = rho * u;
    U(i, 3) = E;

end

%% Movie Setting
% Create the movie handle
F = struct('cdata',[],'colormap',[]);
% Export path for flow field movie
% avi_obj = VideoWriter('test.avi', 'Uncompressed AVI');
if (flag_exp_avi == 1)

    avi_obj = VideoWriter([save_output_folder, 'Sod_Shock_Tube.avi'], 'Motion JPEG AVI');
    avi_obj.Quality = 95;
    avi_obj.FrameRate = 100;
    % open the .avi handle
    open(avi_obj);
    
end

%% Start Calculation & Time marching
% The real cal. time
t = 0;
cnt_step = 0;

while (cnt_step < max_step)
    
    t = t + dt;
    cnt_step = cnt_step + 1;

    % Temporal discretization
    if (flag_tim_mar == 1)
        % 1 - Euler time marching

        if (flag_flu_spl == 1)
            % 1 - Flux Vector Splitting (FVS)
            [F_p, F_n] = Flux_Vect_Split_Common(U, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [xs_new, xt_new, Fh_p, Fh_n, Fx, Fx_p, Fx_n] = Diff_Cons_Common(N, dx, F_p, F_n, flag_spa_typ, flag_upw_typ, flag_scs_typ);
    
        elseif (flag_flu_spl == 2)
            % 2 - Flux Differnce Splitting (FDS)
            [xs_new, xt_new, Fx] = Flux_Diff_Split_Common(U, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

        end

        U = U + (dt * ((-1) * Fx));  % Q(U) = (-1) * Fx

    elseif (flag_tim_mar == 2)
        % 2 - Trapezoid formula (2_od, i.e. improved Euler formula)

        if (flag_flu_spl == 1)
            % 1 - Flux Vector Splitting (FVS)
            [F_p, F_n] = Flux_Vect_Split_Common(U, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx, ~, ~] = Diff_Cons_Common(N, dx, F_p, F_n, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + (dt * ((-1) * Fx));
            [F_p_1, F_n_1] = Flux_Vect_Split_Common(U_1, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [xs_new, xt_new, Fh_p, Fh_n, Fx_1, Fx_p, Fx_n] = Diff_Cons_Common(N, dx, F_p_1, F_n_1, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_1) = (-1) * Fx_1

            U = (0.5 * U) + (0.5 * U_1) + ((0.5 * dt) * ((-1) * Fx_1));
    
        elseif (flag_flu_spl == 2)
            % 2 - Flux Differnce Splitting (FDS)
            [~, ~, Fx] = Flux_Diff_Split_Common(U, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + (dt * ((-1) * Fx));
            [xs_new, xt_new, Fx_1] = Flux_Diff_Split_Common(U_1, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U = (0.5 * U) + (0.5 * U_1) + ((0.5 * dt) * ((-1) * Fx_1));
    
        end

    elseif (flag_tim_mar == 3)
        % 3 - 2_od Runge-Kunta method (Heun formula)

        if (flag_flu_spl == 1)
            % 1 - Flux Vector Splitting (FVS)
            [F_p, F_n] = Flux_Vect_Split_Common(U, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx, ~, ~] = Diff_Cons_Common(N, dx, F_p, F_n, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + ((dt / 3) * ((-1) * Fx));
            [F_p_1, F_n_1] = Flux_Vect_Split_Common(U_1, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx_1, ~, ~] = Diff_Cons_Common(N, dx, F_p_1, F_n_1, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_1) = (-1) * Fx_1

            U_2 = U + (((2 * dt) / 3) * ((-1) * Fx_1));
            [F_p_2, F_n_2] = Flux_Vect_Split_Common(U_2, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [xs_new, xt_new, Fh_p, Fh_n, Fx_2, Fx_p, Fx_n] = Diff_Cons_Common(N, dx, F_p_2, F_n_2, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_2) = (-1) * Fx_2

            U = ((1 / 4) * U) + ((3 / 4) * U_1) + (((3 / 4) * dt) * ((-1) * Fx_2));
    
        elseif (flag_flu_spl == 2)
            % 2 - Flux Differnce Splitting (FDS)
            [~, ~, Fx] = Flux_Diff_Split_Common(U, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + ((dt / 3) * ((-1) * Fx));
            [~, ~, Fx_1] = Flux_Diff_Split_Common(U_1, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_2 = U + (((2 * dt) / 3) * ((-1) * Fx_1));
            [xs_new, xt_new, Fx_2] = Flux_Diff_Split_Common(U_2, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U = ((1 / 4) * U) + ((3 / 4) * U_1) + (((3 / 4) * dt) * ((-1) * Fx_2));
    
        end

    elseif (flag_tim_mar == 4)
        % 4 - 3_od TVD Runge-Kunta method

        if (flag_flu_spl == 1)
            % 1 - Flux Vector Splitting (FVS)
            [F_p, F_n] = Flux_Vect_Split_Common(U, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx, ~, ~] = Diff_Cons_Common(N, dx, F_p, F_n, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + (dt * ((-1) * Fx));
            [F_p_1, F_n_1] = Flux_Vect_Split_Common(U_1, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx_1, ~, ~] = Diff_Cons_Common(N, dx, F_p_1, F_n_1, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_1) = (-1) * Fx_1

            U_2 = ((3 / 4) * U) + ((1 / 4) * U_1) + (((1 * dt) / 4) * ((-1) * Fx_1));
            [F_p_2, F_n_2] = Flux_Vect_Split_Common(U_2, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [xs_new, xt_new, Fh_p, Fh_n, Fx_2, Fx_p, Fx_n] = Diff_Cons_Common(N, dx, F_p_2, F_n_2, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_2) = (-1) * Fx_2

            U = ((1 / 3) * U) + ((2 / 3) * U_2) + (((2 / 3) * dt) * ((-1) * Fx_2));
    
        elseif (flag_flu_spl == 2)
            % 2 - Flux Differnce Splitting (FDS)
            [~, ~, Fx] = Flux_Diff_Split_Common(U, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + (dt * ((-1) * Fx));
            [~, ~, Fx_1] = Flux_Diff_Split_Common(U_1, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_2 = ((3 / 4) * U) + ((1 / 4) * U_1) + (((1 * dt) / 4) * ((-1) * Fx_1));
            [xs_new, xt_new, Fx_2] = Flux_Diff_Split_Common(U_2, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U = ((1 / 3) * U) + ((2 / 3) * U_2) + (((2 / 3) * dt) * ((-1) * Fx_2));
    
        end

    elseif (flag_tim_mar == 5)
        % 5 - 4_od Runge-Kunta method

        if (flag_flu_spl == 1)
            % 1 - Flux Vector Splitting (FVS)
            [F_p, F_n] = Flux_Vect_Split_Common(U, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx, ~, ~] = Diff_Cons_Common(N, dx, F_p, F_n, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + ((1 / 2) * (dt * ((-1) * Fx)));
            [F_p_1, F_n_1] = Flux_Vect_Split_Common(U_1, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx_1, ~, ~] = Diff_Cons_Common(N, dx, F_p_1, F_n_1, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_1) = (-1) * Fx_1

            U_2 = U + ((1 / 2) * (dt * ((-1) * Fx_1)));
            [F_p_2, F_n_2] = Flux_Vect_Split_Common(U_2, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [~, ~, ~, ~, Fx_2, ~, ~] = Diff_Cons_Common(N, dx, F_p_2, F_n_2, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_2) = (-1) * Fx_2

            U_3 = U + (dt * ((-1) * Fx_2));
            [F_p_3, F_n_3] = Flux_Vect_Split_Common(U_3, N, Gamma, Cp, Cv, R, flag_fvs_met);
            [xs_new, xt_new, Fh_p, Fh_n, Fx_3, Fx_p, Fx_n] = Diff_Cons_Common(N, dx, F_p_3, F_n_3, flag_spa_typ, flag_upw_typ, flag_scs_typ); % Q(U_3) = (-1) * Fx_3

            U = ((1 / 3) * (((-1) * U) + U_1 + (2 * U_2) + U_3)) + ((1 / 6) * dt * ((-1) * Fx_3));
    
        elseif (flag_flu_spl == 2)
            % 2 - Flux Differnce Splitting (FDS)
            [~, ~, Fx] = Flux_Diff_Split_Common(U, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_1 = U + ((1 / 2) * (dt * ((-1) * Fx)));
            [~, ~, Fx_1] = Flux_Diff_Split_Common(U_1, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_2 = U + ((1 / 2) * (dt * ((-1) * Fx_1)));
            [~, ~, Fx_2] = Flux_Diff_Split_Common(U_2, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U_3 = U + (dt * ((-1) * Fx_2));
            [xs_new, xt_new, Fx_3] = Flux_Diff_Split_Common(U_3, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ);

            U = ((1 / 3) * (((-1) * U) + U_1 + (2 * U_2) + U_3)) + ((1 / 6) * dt * ((-1) * Fx_3));
    
        end

    end

    if (flag_exp_mov == 1)

    % Save the final solution of U at time t_final
    U_tem = U;
    
    % Cal. the final properties at time t_final
    % Cal. rho, u, E, T, p, c, H according to U
    rho_tem = U_tem(:, 1);
    u_tem = U_tem(:, 2) ./ U_tem(:, 1);
    E_tem = U_tem(:, 3);  % Specific total energy (per unit mass)
    T_tem = ((E_tem ./ rho_tem) - (0.5 .* u_tem .* u_tem)) / Cv;
    p_tem = rho_tem .* R .* T_tem;
    c_tem = sqrt(Gamma * p_tem ./ rho_tem);
    H_tem = (0.5 .* u_tem .* u_tem) + (Cp .* T_tem);  % Specific entropy (per unit mass)
    e_tem = (E_tem ./ rho_tem) - (0.5 .* u_tem .* u_tem);  % Specific internal energy (per unit mass)
    
    % data_tem = analytic_sod(t);  % Call analytic_sod.m function, obtain the analytical solution

    % Monitor properties
        if (cnt_step == 1)
    
            % Plot the basic contour frame (initial fig.)
            h = figure;
            set(0,'CurrentFigure', h)
            set(gcf,'Position',[680,278,660,520]);
            set(gcf,'color','w');
            Plot_Props(t, xp, rho_tem, p_tem, u_tem, e_tem);
    
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            if (flag_exp_mov == 1)
                F = getframe(h);
            end
            
            % Export the movie
            % .avi video
            % Add into the object of pre-created .avi handle
            if (flag_exp_avi == 1)
                writeVideo(avi_obj, F);
            end
            % .gif picture
            % Write in the file with GIF89a format
            if (flag_exp_gif == 1)
                im = frame2im(F);
                [I, map] = rgb2ind(im, 256);
                imwrite(I, map, [save_output_folder, 'Sod_Shock_Tube.gif'], 'GIF', 'Loopcount', inf, 'DelayTime', 0.1);
            end
    
        else
            
            set(0,'CurrentFigure', h)
            Plot_Props(t, xp, rho_tem, p_tem, u_tem, e_tem);
    
            drawnow;
            if (flag_exp_mov == 1)
                F = getframe(h);
            end
    
            % Export the movie
            % .avi video
            % Add into the object of pre-created .avi handle
            if (flag_exp_avi == 1)
                writeVideo(avi_obj, F);
            end
            % .gif picture
            % Write in the file with GIF89a format
            if (flag_exp_gif == 1)
                im = frame2im(F);
                [I, map] = rgb2ind(im, 256);
                imwrite(I, map, [save_output_folder, 'Sod_Shock_Tube.gif'], 'GIF', 'WriteMode', 'append', 'DelayTime', 0.1);
            end
    
        end

    end

end

%% Post-processing

% close the .avi handle
if (flag_exp_avi == 1)
    close(avi_obj);
end

% Save the final solution of U at time t_end
U_end = U;
t_end = t;

% Cal. the final properties at time t_final
% Cal. rho, u, E, T, p, c, H according to U
rho_end = U_end(:, 1);
u_end = U_end(:, 2) ./ U_end(:, 1);
E_end = U_end(:, 3);
T_end = ((E_end ./ rho_end) - (0.5 .* u_end .* u_end)) / Cv;
p_end = rho_end .* R .* T_end;
c_end = sqrt(Gamma * p_end ./ rho_end);
H_end = (0.5 .* u_end .* u_end) + (Cp .* T_end);
e_end = (E_end ./ rho_end) - (0.5 .* u_end .* u_end);

% Call analytic_sod.m function, obtain the analytical solution
data_end = analytic_sod(t_end);

% Visualization
h_end = figure;
% plot(xp, U(:,1)); % rho
set(0,'CurrentFigure', h_end);
set(gcf,'Position',[680,278,660,520]);
set(gcf,'color','w');

% title(zone_title_comb);

subplot(2,2,1);
hold on;
plot(data_end.x, data_end.rho, '-b', 'LineWidth', 1.5);
plot(xp, rho_end, 'bo', 'MarkerSize', 3);
hold off;
legend('Exact', 'Numerical');
xlabel('x');
ylabel('Density');
title(['Plot of Density vs. Position', ' T = ', num2str(t_end, '%.3f'), ' s']);
axis([-0.5, 0.5, 0.0, 1.0]);
grid on;

subplot(2,2,2);
hold on;
plot(data_end.x, data_end.P, '-g', 'LineWidth', 1.5);
plot(xp, p_end, 'go', 'MarkerSize', 3);
hold off;
legend('Exact', 'Numerical');
xlabel('x');
ylabel('Pressure');
title(['Plot of Pressure vs. Position', ' T = ', num2str(t_end, '%.3f'), ' s']);
axis([-0.5, 0.5, 0.0, 1.0]);
grid on;

subplot(2,2,3);
hold on;
plot(data_end.x, data_end.u, '-r', 'LineWidth', 1.5);
plot(xp, u_end, 'ro', 'MarkerSize', 3);
hold off;
legend('Exact', 'Numerical');
xlabel('x');
ylabel('Velocity');
title(['Plot of Velocity vs. Position', ' T = ', num2str(t_end, '%.3f'), ' s']);
axis([-0.5, 0.5, 0.0, 1.0]);
grid on;

subplot(2,2,4);
hold on;
plot(data_end.x, data_end.e, '-m', 'LineWidth', 1.5);
plot(xp, e_end, 'mo', 'MarkerSize', 3);
hold off;
legend('Exact', 'Numerical');
xlabel('x');
ylabel('Specific Internal Energy');
title(['Plot of e vs. Position', ' T = ', num2str(t_end, '%.3f'), ' s']);
axis([-0.5, 0.5, 1.5, 3.0]);
grid on;

% Save the present figure
saveas(gcf, [save_output_folder, 'Results_', zone_title_comb, '_MaxTime_',num2str(max_tot_time, '%.3f'), '.fig']);

% Export the results to Tecplot
% Data preparation
title_cal = zone_title_comb;
zone_title_cal = zone_title_comb;
filename_cal = [save_output_folder, 'Results_Sod_Shock_Tube_', zone_title_cal, '_MaxTime_',num2str(max_tot_time, '%.3f'), '.plt'];
variables_cal = {'X', 'Density', 'Pressure', 'Velocity', 'Specific Internal Energy'};
Mat_Data_cal = [xp(:), rho_end(:), p_end(:), u_end(:), e_end(:)];
IJK_cal = length(xp);

% Create the file
if exist(filename_cal, 'file') 
    delete(filename_cal)
end
f_id_cal = fopen(filename_cal, 'a');
fclose(f_id_cal);

% Export the exact solution to Tecplot
% Create the header
plt_Head(filename_cal, title_cal, variables_cal);
% Create the format of zone(point)
plt_Zone(filename_cal, zone_title_cal, IJK_cal, Mat_Data_cal);

title_ana = 'Exact';
zone_title_ana = 'Exact';
filename_ana = [save_output_folder, 'Results_Sod_Shock_Tube_', zone_title_ana, '_MaxTime_',num2str(max_tot_time, '%.3f'), '.plt'];
variables_ana = {'X', 'Density', 'Pressure', 'Velocity', 'Specific Internal Energy'};
Mat_Data_ana = [data_end.x(:), data_end.rho(:), data_end.P(:), data_end.u(:), data_end.e(:)];
IJK_ana = length(data_end.x);

% Create the file
if exist(filename_ana, 'file') 
    delete(filename_ana)
end
f_id_ana = fopen(filename_ana, 'a');
fclose(f_id_ana);

% Create the header
plt_Head(filename_ana, title_ana, variables_ana);
% Create the format of zone(point)
plt_Zone(filename_ana, zone_title_ana, IJK_ana, Mat_Data_ana);

% Save data of the program
save([save_output_folder, 'Results_Variables_', zone_title_comb, '_MaxTime_',num2str(max_tot_time, '%.3f'), '.mat']);

% END OF THE PROGRAM!

    % Function module: Create the header
    function plt_Head(filename, title, variables)
    
        % Create the header
        f_id = fopen(filename, 'a');
        % Name
        if ~isempty(title)
            s = ['TITLE = "', title, '"'];
            fprintf(f_id, '%s \r\n', s);
        end
        % Variables
        v = numel(variables);
        s = 'VARIABLES = ';
        for k = 1 : v
            if k ~= 1
                s = [s, ','];
            end
            s = [s, ' "', variables{k}, '"'];
        end
        fprintf(f_id, '%s \r\n', s);
    
        fclose(f_id);
    
    end
    
    % Function module: Create the format of zone(point)
    function plt_Zone(filename, zone_title, IJK, Mat_Data)
        
        % Create the format of zone(point)
        f_id = fopen(filename, 'a');
        N = size(Mat_Data, 1);
        
        Dim = numel(IJK);
        if (Dim == 1)
            s = ['zone I =', num2str(IJK(1))];
        elseif (Dim == 2)
            s = ['zone I =', num2str(IJK(1)), ', J =', num2str(IJK(2))];
        elseif (Dim == 3)
            s = ['zone I =', num2str(IJK(1)), ', J =', num2str(IJK(2)), ', K =', num2str(IJK(3))];
        end
        
        % Title
        if ~isempty(zone_title)
            s = [s, ', t = "', zone_title, '"'];
        end
        fprintf(f_id, '%s \r\n', s);
        % Point format
        s = 'DATAPACKING = point';
        fprintf(f_id, '%s \r\n', s);
        
        % Data introduction
        for k = 1 : N
            fprintf(f_id, '%s \r\n', num2str(Mat_Data(k,:)));
        end
        
        fclose(f_id);
    
    end

