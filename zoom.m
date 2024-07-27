function [alpha_star, f_eval, g_eval] = zoom(alpha_low, alpha_high, Phi, control)
    Phi_0 = double(subs(Phi, 0));  % Phi(0)
    Phi_prime = gradient(Phi);  % Derivative of Phi
    Phi_prime_0 = double(subs(Phi_prime, 0));  % Derivative of Phi at alpha = 0
    f_eval = 1;  % Function evaluations
    g_eval = 1;  % Gradient evaluations
    c1 = 1e-4;
    if control == 1
        c2 = 0.9;  % For Newton
    else
        c2 = 0.1;  % For steepest descent
    end
    flag = 1;
    while flag
        alpha = 0.5*(alpha_low + alpha_high);
        Phi_i = double(subs(Phi, alpha));
        f_eval = f_eval + 1;
        if (sign(alpha_low-alpha_high))*double(subs(Phi_prime, alpha_low)) <= 1e-5
            g_eval = g_eval + 1;
            alpha_star = alpha;
            flag = 0;
            break;
        end
        if Phi_i > (Phi_0 + c1*alpha*Phi_prime_0) || Phi_i >= double(subs(Phi, alpha_low))
            f_eval = f_eval + 1;
            alpha_high = alpha;
        else
            Phi_prime_i = double(subs(Phi_prime, alpha));
            g_eval = g_eval + 1;
            if abs(Phi_prime_i) <= (-c2 *Phi_prime_0)
                alpha_star = alpha;
                flag = 0;
                break;
            end
            if Phi_prime_i*(alpha_high-alpha_low) >= 0
                alpha_high = alpha_low;
            end
            alpha_low = alpha;
        end
    end
end