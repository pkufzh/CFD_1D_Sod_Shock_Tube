
% Function Module: Flux Vector Splitting (FVS) with different methods

function [F_p, F_n] = Flux_Vect_Split_Common(U, N, Gamma, Cp, Cv, R, flag_fvs_met)

    if (flag_fvs_met == 1)
        % 1 - FVS - Steger-Warming (S-W)

        % Initiation
        lambda = zeros(3, 1);    % eigenvalues matrix 
        lambda_p = zeros(3, 1);  % pos. eigenvalues matrix
        lambda_n = zeros(3, 1);  % neg. eigenvalues matrix
    
        F_p = zeros(N, 3);   % pos. invisid vector
        F_n = zeros(N, 3);   % neg. invisid vector
    
        em = 1e-3;

        % Splitting
        for i = 1 : N
    
            % Step 1: Cal. rho, u, p, T, c according to U
            rho = U(i, 1);
            u = U(i, 2) / U(i, 1);
            E = U(i, 3);
            T = ((E / rho) - (0.5 * u * u)) / Cv;
            p = rho * R * T;
            c = sqrt(Gamma * p / rho);
            
            % Step 2: Cal. eigenvalues lamba
            lambda(1) = u;
            lambda(2) = u - c;
            lambda(3) = u + c;
    
            % Step 3: Splitting the eigenvalues lamba into lamba (S-W)
            for k = 1 : 3
                lambda_p(k) = (lambda(k) + sqrt((lambda(k)^2) + (em^2))) / 2;
                lambda_n(k) = (lambda(k) - sqrt((lambda(k)^2) + (em^2))) / 2;
            end
    
            % Step 4: Cal. F+ & F-
            w_p = ((3 - Gamma) * (lambda_p(2) + lambda_p(3)) * c * c) / (2 * (Gamma - 1));
            F_p(i, 1) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_p(1)) + lambda_p(2) + lambda_p(3));
            F_p(i, 2) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_p(1) * u) + (lambda_p(2) * (u - c)) + (lambda_p(3) * (u + c)));
            % !!!
            F_p(i, 3) = (rho / (2 * Gamma)) * (((Gamma - 1) * lambda_p(1) * u * u) + ((lambda_p(2) / 2) * (u - c) * (u - c)) + (((lambda_p(3) / 2) * (u + c) * (u + c)) + w_p));
    
            w_n = ((3 - Gamma) * (lambda_n(2) + lambda_n(3)) * c * c) / (2 * (Gamma - 1));
            F_n(i, 1) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_n(1)) + lambda_n(2) + lambda_n(3));
            F_n(i, 2) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_n(1) * u) + (lambda_n(2) * (u - c)) + (lambda_n(3) * (u + c)));
            % !!!
            F_n(i, 3) = (rho / (2 * Gamma)) * (((Gamma - 1) * lambda_n(1) * u * u) + ((lambda_n(2) / 2) * (u - c) * (u - c)) + (((lambda_n(3) / 2) * (u + c) * (u + c)) + w_n));
    
        end

    elseif (flag_fvs_met == 2)
        % 2 - FVS - Lax-Friedrich (L-F)

        % Initiation
        lambda = zeros(3, 1);     % eigenvalues matrix 
        lambda_p = zeros(3, 1);  % pos. eigenvalues matrix
        lambda_n = zeros(3, 1);  % neg. eigenvalues matrix
    
        F_p = zeros(N, 3);   % pos. invisid vector
        F_n = zeros(N, 3);   % neg. invisid vector
        
        % Splitting

        lambda_s_global = 0;
        for i = 1 : N

            rho = U(i, 1);
            u = U(i, 2) / U(i, 1);
            E = U(i, 3);
            T = ((E / rho) - (0.5 * u * u)) / Cv;
            p = rho * R * T;
            c = sqrt(Gamma * p / rho);
            
            % Global splitting - plus eigenvalue
            lambda_s_global = max(lambda_s_global, abs(u) + c);

        end

        for i = 1 : N
    
            % Step 1: Cal. rho, u, p, T, c according to U
            rho = U(i, 1);
            u = U(i, 2) / U(i, 1);
            E = U(i, 3);
            T = ((E / rho) - (0.5 * u * u)) / Cv;
            p = rho * R * T;
            c = sqrt(Gamma * p / rho);

            % Step 2: Cal. eigenvalues lamba
            lambda(1) = u;
            lambda(2) = u - c;
            lambda(3) = u + c;

            lambda_s_local = abs(u) + c;

            % Step 3: Splitting the eigenvalues lamba into lamba (L-F)
            % Local splitting - plus eigenvalue or Global
            for k = 1 : 3
                lambda_p(k) = (lambda(k) + lambda_s_local) / 2;
                lambda_n(k) = (lambda(k) - lambda_s_local) / 2;
            end
    
            % Step 4: Cal. F+ & F-
            w_p = ((3 - Gamma) * (lambda_p(2) + lambda_p(3)) * c * c) / (2 * (Gamma - 1));
            F_p(i, 1) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_p(1)) + lambda_p(2) + lambda_p(3));
            F_p(i, 2) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_p(1) * u) + (lambda_p(2) * (u - c)) + (lambda_p(3) * (u + c)));
            F_p(i, 3) = (rho / (2 * Gamma)) * (((Gamma - 1) * lambda_p(1) * u * u) + ((lambda_p(2) / 2) * (u - c) * (u - c)) + (((lambda_p(3) / 2) * (u + c) * (u + c)) + w_p));
    
            w_n = ((3 - Gamma) * (lambda_n(2) + lambda_n(3)) * c * c) / (2 * (Gamma - 1));
            F_n(i, 1) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_n(1)) + lambda_n(2) + lambda_n(3));
            F_n(i, 2) = (rho / (2 * Gamma)) * ((2 * (Gamma - 1) * lambda_n(1) * u) + (lambda_n(2) * (u - c)) + (lambda_n(3) * (u + c)));
            F_n(i, 3) = (rho / (2 * Gamma)) * (((Gamma - 1) * lambda_n(1) * u * u) + ((lambda_n(2) / 2) * (u - c) * (u - c)) + (((lambda_n(3) / 2) * (u + c) * (u + c)) + w_n));

        end

    elseif (flag_fvs_met == 3)
        % 3 - Van Leer

        % Initiation
        F = zeros(1, 3);     % invisid vector
        F_p = zeros(N, 3);   % pos. invisid vector
        F_n = zeros(N, 3);   % neg. invisid vector

        % Splitting
        for i = 1 : N
    
            % Step 1: Cal. rho, u, p, T, c, F according to U
            rho = U(i, 1);
            u = U(i, 2) / U(i, 1);
            E = U(i, 3);
            T = ((E / rho) - (0.5 * u * u)) / Cv;
            p = rho * R * T;
            c = sqrt(Gamma * p / rho);

            % F(1) = rho * u;
            % F(2) = (rho * u * u) + p;
            % F(3) = u * (E + p);

            F(1) = U(i, 2);
            F(2) = ((Gamma - 1) * U(i, 3)) + (((3 - Gamma) / 2) * ((U(i, 2) * U(i, 2)) / U(i, 1)));
            F(3) = (Gamma * (U(i, 2) * U(i, 3)) / U(i, 1)) + (((Gamma - 1) / 2) * (U(i, 2)^3 / U(i, 1)^2));


            % Step 2: Cal. Mach number Ma
            Ma = u / c;

            % Step 3: Splitting by discuss different Ma cases
            if (Ma >= 1)
                F_p(i, :) = F(1, :);
                F_n(i, :) = zeros(1, 3);
            elseif (Ma <= -1)
                F_p(i, :) = zeros(1, 3);
                F_n(i, :) = F(1, :);
            else
                F1_p = rho * c * (((Ma + 1) / 2)^2);
                F1_n = (-1) * rho * c * (((Ma - 1) / 2)^2);

                F_p(i, 1) = F1_p;
                F_p(i, 2) = (F1_p / Gamma) * (((Gamma - 1) * u) + (2 * c));
                F_p(i, 3) = (F1_p / (2 * (Gamma^2 - 1))) * (((Gamma - 1) * u) + (2 * c))^2;

                F_n(i, 1) = F1_n;
                F_n(i, 2) = (F1_n / Gamma) * (((Gamma - 1) * u) - (2 * c));
                F_n(i, 3) = (F1_n / (2 * (Gamma^2 - 1))) * (((Gamma - 1) * u) - (2 * c))^2;
            end

        end

    elseif (flag_fvs_met == 4)
        % 4 - Liou-Steffen Splitting - Advection Upstream Splitting (AUSM) Method

        % Initiation
        Fc_p = zeros(1, 3);  % pos. invisid vector convective part
        Fc_n = zeros(1, 3);  % neg. invisid vector convective part
        Fp_p = zeros(1, 3);  % pos. invisid vector pressure part
        Fp_n = zeros(1, 3);  % neg. invisid vector pressure part
        F_p = zeros(N, 3);   % pos. invisid vector
        F_n = zeros(N, 3);   % neg. invisid vector

        % Splitting
        for i = 1 : N
    
            % Step 1: Cal. rho, u, p, T, c, F according to U
            rho = U(i, 1);
            u = U(i, 2) / U(i, 1);
            E = U(i, 3);
            T = ((E / rho) - (0.5 * u * u)) / Cv;
            p = rho * R * T;
            c = sqrt(Gamma * p / rho);
            H = (0.5 * u * u) + (Cp * T);

            % Step 2: Cal. Mach number Ma
            Ma = u / c;

            % Step 3: Splitting by discuss different Ma cases
            if (Ma > 1)
                Ma_p = Ma;
                Ma_n = 0;
                p_p = p;
                p_n = 0;
            elseif (Ma < -1)
                Ma_p = 0;
                Ma_n = Ma;
                p_p = 0;
                p_n = p;
            else
                Ma_p = ((Ma + 1)^2) / 4;
                Ma_n = ((-1) * ((Ma - 1)^2)) / 4;
                p_p = p * ((1 + Ma) / 2);
                p_n = p * ((1 - Ma) / 2);
            end

            Fc_p(1) = rho * c * Ma_p;
            Fc_p(2) = rho * c * Ma_p * u;
            Fc_p(3) = rho * c * Ma_p * H;

            Fc_n(1) = rho * c * Ma_n;
            Fc_n(2) = rho * c * Ma_n * u;
            Fc_n(3) = rho * c * Ma_n * H;

            Fp_p(1) = 0;
            Fp_p(2) = p_p;
            Fp_p(3) = 0;

            Fp_n(1) = 0;
            Fp_n(2) = p_n;
            Fp_n(3) = 0;

            F_p(i, 1) = Fc_p(1) + Fp_p(1);
            F_p(i, 2) = Fc_p(2) + Fp_p(2);
            F_p(i, 3) = Fc_p(3) + Fp_p(3);

            F_n(i, 1) = Fc_n(1) + Fp_n(1);
            F_n(i, 2) = Fc_n(2) + Fp_n(2);
            F_n(i, 3) = Fc_n(3) + Fp_n(3);

        end

    end

end

