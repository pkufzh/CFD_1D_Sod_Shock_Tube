
% Function Module: the approach to cal. difference for Fx from F_p and F_n with conservation form
% Note: the upwind schemes are converted into conservative form

function [xs_new, xt_new, Fh_p, Fh_n, Fx, Fx_p, Fx_n] = Diff_Cons_Common(N, dx, F_p, F_n, flag_spa_typ, flag_upw_typ, flag_scs_typ)

    Fx = zeros(N, 3);    % invisid term
    Fx_p = zeros(N, 3);  % pos. invisid term
    Fx_n = zeros(N, 3);  % neg. invisid term

    Fh_p = zeros(N, 3); % half point (j + 1/2) pos. invisid vector
    Fh_n = zeros(N, 3); % half point (j + 1/2) neg. invisid vector
    Fh = zeros(N, 3); % half point (j + 1/2) invisid vector

    % Step 5: Cal. the flux derivative with different schemes (& shock-capturing)
    if (flag_spa_typ == 1)
        % 1 - General Upwind & Compact Schemes (forward / backward)

        if (flag_upw_typ == 1)
            % 1 - 1_od upwind scheme (2 points)
            ks = -1;           % the start corner index relative to j
            kt = 0;
            kn = kt - ks + 1;  % number of the coefficients
            kp = (ks : kt);
            
            a = [];
            a(1) = -1;
            a(2) = 1;
            
            b = zeros((kn - 1), 1);
            b(1) = (-1) * a(1);
            for k = 2 : (kn - 1)
                b(k) = b(k - 1) - a(k);
            end

        elseif (flag_upw_typ == 2)
            % 2 - 2_od upwind scheme (3 points)
            ks = -2;
            kt = 0;
            kn = kt - ks + 1;
            kp = (ks : kt);
            
            a = [];
            a(1) = 1 / 2;
            a(2) = -4 / 2;
            a(3) = 3 / 2;
            
            b = zeros((kn - 1), 1);
            b(1) = (-1) * a(1);
            for k = 2 : (kn - 1)
                b(k) = b(k - 1) - a(k);
            end

        elseif (flag_upw_typ == 3)
            % 3 - 3_od upwind scheme (4 points with bias)
            ks = -2;
            kt = 1;
            kn = kt - ks + 1;  % number of the coefficients
            kp = (ks : kt);
            
            a = [];
            a(1) = 1 / 6;
            a(2) = -6 / 6;
            a(3) = 3 / 6;
            a(4) = 2 / 6;
            
            b = zeros((kn - 1), 1);
            b(1) = (-1) * a(1);
            for k = 2 : (kn - 1)
                b(k) = b(k - 1) - a(k);
            end

        elseif (flag_upw_typ == 4)
            % 4 - 5_od upwind scheme (6 points with bias)
            ks = -3;
            kt = 2;
            kn = kt - ks + 1;  % number of the coefficients
            kp = (ks : kt);
            
            a = [];
            a(1) = -2 / 60;
            a(2) = 15 / 60;
            a(3) = -60 / 60;
            a(4) = 20 / 60;
            a(5) = 30 / 60;
            a(6) = -3 / 60;

            b = zeros((kn - 1), 1);
            b(1) = (-1) * a(1);
            for k = 2 : (kn - 1)
                b(k) = b(k - 1) - a(k);
            end

        end
        
        % [Core algorithm]
        % Uniform cal. procedure according to the coeff.
        xs = 1 + abs(kp(2));
        xt = N - abs(kp(2));
        for j = xs : xt
            for k = 1 : (kn - 1)
                Fh_p(j, :) = Fh_p(j, :) + b(k) * F_p(j + kp(k + 1), :);  % F+ (j + 1/2)
                Fh_n(j, :) = Fh_n(j, :) + b(k) * F_n(j - kp(k + 1), :);  % F- (j - 1/2)  different points?
            end
        end

        xs_new = xs + 1;
        xt_new = xt - 1;
        for j = xs_new : xt_new
            Fx_p(j, :) = (Fh_p(j, :) - Fh_p(j - 1, :)) / dx;
            Fx_n(j, :) = (Fh_n(j + 1, :) - Fh_n(j, :)) / dx;
            Fx(j, :) = Fx_p(j, :) + Fx_n(j, :);
        end

    elseif (flag_spa_typ == 2)
        % 2 - Special Shock-Capturing Schemes

        if (flag_scs_typ == 1)
            % 1 - TVD - Total Variation Diminishing Scheme (Van Leer Limiter)
            xs = 2;
            xt = N - xs;

            em = 1e-5;

            for j = xs : xt
                % Using Van Leer limiter (or other limiters)
                r_p = (F_p(j, :) - F_p(j - 1, :)) ./ (F_p(j + 1, :) - F_p(j, :) + em);    % NaN? em!
                r_n = (F_n(j + 2, :) - F_n(j + 1, :)) ./ (F_n(j + 1, :) - F_n(j, :) + em);
                Phi_p = (r_p + abs(r_p)) ./ (1 + r_p);
                Phi_n = (r_n + abs(r_n)) ./ (1 + r_n);

                Fh_p(j, :) = F_p(j, :) + 0.5 * (Phi_p .* (F_p(j + 1, :) - F_p(j, :)));
                Fh_n(j, :) = F_n(j + 1, :) - 0.5 * (Phi_n .* (F_n(j + 1, :) - F_n(j, :)));
            end
            
            % start point revision?
            xs_new = xs + 1;
            xt_new = xt;

            for j = xs_new : xt_new
                Fx_p(j, :) = (Fh_p(j, :) - Fh_p(j - 1, :)) / dx;
                Fx_n(j, :) = (Fh_n(j, :) - Fh_n(j - 1, :)) / dx;
                Fx(j, :) = Fx_p(j, :) + Fx_n(j, :);
            end
            
        elseif (flag_scs_typ == 2)
            % 2 - NND - Non-oscillatory, Non-free-paremeters Dissipative Difference Scheme
            xs = 2;
            xt = N - xs;

            for j = xs : xt
                Fh_p(j, :) = F_p(j, :) +     (0.5 * Cal_Minmod((F_p(j, :) - F_p(j - 1, :)), (F_p(j + 1, :) - F_p(j, :))));
                Fh_n(j, :) = F_n(j + 1, :) - (0.5 * Cal_Minmod((F_n(j + 1, :) - F_n(j, :)), (F_n(j + 2, :) - F_n(j + 1, :))));
                Fh(j, :) = Fh_p(j, :) + Fh_n(j, :);
            end

            xs_new = xs + 1;
            xt_new = xt;

            for j = xs_new : xt_new
                Fx(j, :) = (Fh(j, :) - Fh(j - 1, :)) / dx;
            end

        elseif (flag_scs_typ == 3)
            % 3 - WENO - Weighted Essentially Non-Oscillatory Method
            % (Jiang & Shu, 1996) 5 order WENO scheme

            % Parameters
            C = zeros(1, 3);

            C(1) = 1 / 10;
            C(2) = 6 / 10;
            C(3) = 3 / 10;

            p = 2;
            em = 1e-6;

            % a > 0 case     
            beta_p = zeros(3, 3);
            alpha_p = zeros(3, 3);
            omega_p = zeros(3, 3); % Weights
            Fh_p_c = zeros(3, 3);

            % Fh_p_1 = zeros(N, 3);  % Stencil 1
            % Fh_p_2 = zeros(N, 3);  % Stencil 2
            % Fh_p_3 = zeros(N, 3);  % Stencil 3
            Fh_p = zeros(N, 3);  % Sum of weighted stencils (number = 3)

            % a < 0 case
            beta_n = zeros(3, 3);
            alpha_n = zeros(3, 3);
            omega_n = zeros(3, 3); % Weights
            Fh_n_c = zeros(3, 3);
            Fh_n = zeros(N, 3);  % Sum of weighted stencils (number = 3)

            xs = 3;
            xt = N - xs + 1;
            
            for j = xs : xt

                % a > 0, pos. flux case
                % Cal. weights of each stencil !!!
                beta_p(1, :) = ((1 / 4) * (F_p(j - 2, :) - 4 * F_p(j - 1, :) + 3 * F_p(j, :)).^2) + ((13 / 12) * (F_p(j - 2, :) - 2 * F_p(j - 1, :) + F_p(j, :)).^2);
                beta_p(2, :) = ((1 / 4) * (F_p(j - 1, :) - F_p(j + 1, :)).^2) + ((13 / 12) * (F_p(j - 1, :) - 2 * F_p(j, :) + F_p(j + 1, :)).^2);
                beta_p(3, :) = ((1 / 4) * (3 * F_p(j, :) - 4 * F_p(j + 1, :) + F_p(j + 2, :)).^2) + ((13 / 12) * (F_p(j, :) - 2 * F_p(j + 1, :) + F_p(j + 2, :)).^2);
                
                for k = 1 : 3
                    alpha_p(k, :) = C(k) ./ ((em + beta_p(k, :)).^p);
                end
                alpha_p_sum = sum(alpha_p);

                for k = 1 : 3
                    omega_p(:, k) = alpha_p(:, k) / alpha_p_sum(k); % !!! notice orders
                end
                
                % Construct stencils !!! (F+ (j + 1/2))
                Fh_p_c(1, :) = ((1 / 3) * F_p(j - 2, :)) - ((7 / 6) * F_p(j - 1, :)) + ((11 / 6) * F_p(j, :));
                Fh_p_c(2, :) = ((-1) * (1 / 6) * F_p(j - 1, :)) + ((5 / 6) * F_p(j, :)) + ((1 / 3) * F_p(j + 1, :));
                Fh_p_c(3, :) = ((1 / 3) * F_p(j, :)) + ((5 / 6) * F_p(j + 1, :)) - ((1 / 6) * F_p(j + 2, :));

                % omega_p' * Fh_p_c
                omega_p_c = omega_p';

                % Sum up all the stencils with weight to obtain F+ (j + 1/2)
                for k = 1 : 3
                    Fh_p(j, k) = (omega_p_c(k, :) * Fh_p_c(:, k));
                end

                % a < 0, neg. flux case
                % Cal. weights of each stencil !!!
                beta_n(1, :) = ((1 / 4) * (F_n(j + 2, :) - 4 * F_n(j + 1, :) + 3 * F_n(j, :)).^2) + ((13 / 12) * (F_n(j + 2, :) - 2 * F_n(j + 1, :) + F_n(j, :)).^2);
                beta_n(2, :) = ((1 / 4) * (F_n(j + 1, :) - F_n(j - 1, :)).^2) + ((13 / 12) * (F_n(j + 1, :) - 2 * F_n(j, :) + F_n(j - 1, :)).^2);
                beta_n(3, :) = ((1 / 4) * (3 * F_n(j, :) - 4 * F_n(j - 1, :) + F_n(j - 2, :)).^2) + ((13 / 12) * (F_n(j, :) - 2 * F_n(j - 1, :) + F_n(j - 2, :)).^2);
                
                for k = 1 : 3
                    alpha_n(k, :) = C(k) ./ ((em + beta_n(k, :)).^p);
                end
                alpha_n_sum = sum(alpha_n);

                for k = 1 : 3
                    omega_n(:, k) = alpha_n(:, k) / alpha_n_sum(k);
                end
                
                % Construct stencils !!! (F- (j - 1/2))
                Fh_n_c(1, :) = ((1 / 3) * F_n(j + 2, :)) - ((7 / 6) * F_n(j + 1, :)) + ((11 / 6) * F_n(j, :));
                Fh_n_c(2, :) = ((-1) * (1 / 6) * F_n(j + 1, :)) + ((5 / 6) * F_n(j, :)) + ((1 / 3) * F_n(j - 1, :));
                Fh_n_c(3, :) = ((1 / 3) * F_n(j, :)) + ((5 / 6) * F_n(j - 1, :)) - ((1 / 6) * F_n(j - 2, :));

                % omega_p' * Fh_p_c
                omega_n_c = omega_n';

                % Sum up all the stencils with weight to obtain F- (j - 1/2)
                for k = 1 : 3
                    Fh_n(j, k) = (omega_n_c(k, :) * Fh_n_c(:, k));
                end

            end

            xs_new = xs + 1;
            xt_new = xt - 1;

            for j = xs_new : xt_new
                Fx_p(j, :) = (Fh_p(j, :) - Fh_p(j - 1, :)) / dx;
                Fx_n(j, :) = (Fh_n(j + 1, :) - Fh_n(j, :)) / dx;
                Fx(j, :) = Fx_p(j, :) + Fx_n(j, :);
            end            

        end
        
    end

end

