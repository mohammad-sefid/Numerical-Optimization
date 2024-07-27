function [alpha_star, f_eval, g_eval] = linesearch(Phi, control)
    alpha_prev = 0;
    alpha_max = 5;
    alpha = 0.5*(alpha_prev + alpha_max);
    c1 = 1e-4;
    if control == 1
        c2 = 0.9;  % For Newton
    else
        c2 = 0.1;  % For steepest descent
    end
    Phi_0 = double(subs(Phi, 0));  % Phi(0)
    Phi_prime = gradient(Phi);  % Derivative of Phi
    Phi_prime_0 = double(subs(Phi_prime, 0));  % Derivative of Phi at alpha = 0
    flag = 1;
    f_eval = 1;  % Function evaluations
    g_eval = 1;  % Gradient evaluations
    while flag
        Phi_i = double(subs(Phi, alpha));
        Phi_prev = double(subs(Phi, alpha_prev));
        f_eval = f_eval + 2;
        if Phi_i > (Phi_0 + c1*alpha*Phi_prime_0) || (Phi_i > Phi_prev && f_eval > 2)
            [alpha_star, f_e, g_e] = zoom(alpha_prev, alpha, Phi, control);
            f_eval = f_eval + f_e;
            g_eval = g_eval + g_e;
            flag = 0;
            break;
        else
            Phi_prime_i = double(subs(Phi_prime, alpha));
            g_eval = g_eval + 1;
            if abs(Phi_prime_i) <= (-c2 * Phi_prime_0)
                alpha_star = alpha;
                flag = 0;
                break;
            end
            if Phi_prime_i >= 0
                [alpha_star, f_e, g_e] = zoom(alpha, alpha_prev, Phi, control);
                f_eval = f_eval + f_e;
                g_eval = g_eval + g_e;
                flag = 0;
                break;
            else
                alpha_prev = alpha;
                alpha = 0.5*(alpha + alpha_max);
            end
        end
    end
end