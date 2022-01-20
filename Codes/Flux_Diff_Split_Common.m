
% Function Module: Flux Difference Splitting (FDS) with different methods

function [xs_new, xt_new, Fx] = Flux_Diff_Split_Common(U, N, dx, Gamma, Cp, Cv, R, flag_fds_met, flag_spa_typ, flag_upw_typ, flag_scs_typ)

    if (flag_fds_met == 1)
        % FDS - Roe scheme
        % The variables of all grid points at time n are known.
        Fh_l = zeros(N, 3);
        Fh_r = zeros(N, 3);
        U_ave = zeros(N, 3);  % Roe ave. bar(U) martix
        F_ave = zeros(N, 3);  % Roe ave. F(bar(U)) martix
        Fh = zeros(N, 3);
        Fx = zeros(N, 3);
        A_ave = zeros(3, 3);  % Roe ave. Jacobian A(bar(U)) martix

        em = 1e-5;

        % Step 1：Cal. Ur, Ul with difference schemes
        % Cal. all the variables of all grid points
        [xs, xt, Uh_l, Uh_r, ~, ~, ~] = Diff_Cons_Common(N, dx, U, U, flag_spa_typ, flag_upw_typ, flag_scs_typ);

        % !!! Special case: handling Upwind schemes & WENO schemes (F_l(j + 1/2), F_n(j - 1/2))
        if (flag_spa_typ == 1) || ((flag_spa_typ == 2) && (flag_scs_typ == 3))

            xt = xt - 1;
            % Shifting Uh_r from (j - 1/2) to (j + 1/2)
            for j = xs : xt
                Uh_r(j, :) = Uh_r(j + 1, :);
            end
            
        end

        for j = xs : xt

            % Step 2: Cal. the Roe average value bar(U) with Roe formula
            rho_l = Uh_l(j, 1);
            u_l = Uh_l(j, 2) / Uh_l(j, 1);
            E_l = Uh_l(j, 3);
            T_l = ((E_l / rho_l) - (0.5 * u_l * u_l)) / Cv;
            p_l = rho_l * R * T_l;
            % c_l = sqrt(Gamma * p_l / rho_l);
            H_l = (0.5 * u_l * u_l) + (Cp * T_l);  % Or H = (E + p) / rho
            
            % Cal. F(Ul)
            Fh_l(j, 1) = rho_l * u_l;
            Fh_l(j, 2) = (rho_l * u_l * u_l) + p_l;
            Fh_l(j, 3) = u_l * (E_l + p_l);

            rho_r = Uh_r(j, 1);
            u_r = Uh_r(j, 2) / Uh_r(j, 1);
            E_r = Uh_r(j, 3);
            T_r = ((E_r / rho_r) - (0.5 * u_r * u_r)) / Cv;
            p_r = rho_r * R * T_r;
            % c_r = sqrt(Gamma * p_r / rho_r);
            H_r = (0.5 * u_r * u_r) + (Cp * T_r);
            
            % Cal. F(Ur)
            Fh_r(j, 1) = rho_r * u_r;
            Fh_r(j, 2) = (rho_r * u_r * u_r) + p_r;
            Fh_r(j, 3) = u_r * (E_r + p_r);

            % Cal. ave. bar(U) martix
            rho_ave = ((sqrt(rho_l) + sqrt(rho_r)) / 2)^2;
            u_ave = ((sqrt(rho_l) * u_l + sqrt(rho_r) * u_r) / (2 * sqrt(rho_ave)));
            H_ave = ((sqrt(rho_l) * H_l + sqrt(rho_r) * H_r) / (2 * sqrt(rho_ave)));
            p_ave = ((Gamma - 1) / Gamma) * ((rho_ave * H_ave) - (0.5 * rho_ave * u_ave * u_ave));
            c_ave = sqrt((Gamma - 1) * (H_ave - (0.5 * u_ave * u_ave)));
            E_ave = (rho_ave * H_ave) - p_ave;

            U_ave(j, 1) = rho_ave;
            U_ave(j, 2) = rho_ave * u_ave;
            U_ave(j, 3) = E_ave;

            % Cal. ave. F(bar(U)) martix
            F_ave(j, 1) = rho_ave * u_ave;
            F_ave(j, 2) = (rho_ave * u_ave * u_ave) + p_ave;
            F_ave(j, 3) = u_ave * (E_ave + p_ave);

            % Cal. the Jacobian matrix A(bar(U)) using the relation F(U) = AU (Invalid! Jacobian quality fails!)
            % Or using the direct cal. (Valid)
            % A_ave = (F_ave(j, :)') / (U_ave(j, :)');  % Martix cal. failed

            A_ave(1, 1) = 0;
            A_ave(1, 2) = 1;
            A_ave(1, 3) = 0;
            A_ave(2, 1) = (-1) * ((3 - Gamma) / 2) * u_ave * u_ave;
            A_ave(2, 2) = (3 - Gamma) * u_ave;
            A_ave(2, 3) = Gamma - 1;
            A_ave(3, 1) = (((Gamma - 2) / 2) * u_ave * u_ave * u_ave) - ((u_ave * c_ave * c_ave) / (Gamma - 1));
            A_ave(3, 2) = ((c_ave * c_ave) / (Gamma - 1)) + (((3 - Gamma) / 2) * u_ave * u_ave);
            A_ave(3, 3) = Gamma * u_ave;

            % Step 3: Cal. the Jacobian matrix A(bar(U))
            [V, G] = eig(A_ave);
            S = inv(V);
            % A_ave = inv(S) * G * S;  % Confirm

            % Step 4: Cal. the absolute Jacobian matrix abs(A(bar(U)))
            G_abs = zeros(3, 3);
            % Entropy correction
            for i = 1 : 3
                if (abs(G(i, i)) > em)
                    G_abs(i, i) = abs(G(i, i));
                else
                    G_abs(i, i) = ((G(i, i) * G(i, i)) + (em * em)) / (2 * em);
                end
            end
            A_ave_abs = inv(S) * G_abs * S;

            % Step 5: Cal. F(j + 1/2)
            Fh(j, :) = ((0.5 * (Fh_r(j, :)' + Fh_l(j, :)')) - (0.5 * A_ave_abs * (Uh_r(j, :)' - Uh_l(j, :)')))';

        end

        xs_new = xs + 1;
        xt_new = xt;

        for j = xs_new : xt_new

            % Step 6: Cal. the spatial derivative Fx_j
            Fx(j, :) = (Fh(j, :) - Fh(j - 1, :)) / dx;

        end

    end

end

