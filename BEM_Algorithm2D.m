%% BEM Method
% To be done for each element/radius 

function [pn_r,pt_r, alpha_mem]= BEM_Algorithm2D(V0, R, B, rho, n_airfoil, lambda, theta_p, r, c, beta, k, aoa_list, Cl_list, Cd_list, Thick_prof_list, itt)
    
    % BEM calculation data
    eps = 1e-5;                % [-] : Convergence criterium
    BetzLimit = 1/3;            % [-] : Betz limit

    % Initialising a (=an) and a' (=at)
    n = 1;                       % [-] : Iteration step number
    an_old = 0;             % [-] : (Normal) Induction factor at step n-1
    at_old = 0;             % [-] : Tangential induction factor at step n-1
    an_new = 0;                 % [-] : (Normal) Induction factor at step n
    at_new = 0;                 % [-] : Tangential induction factor at step n

    alpha_mem = 0;

    running = true;

    while running
        % 
        an_old = an_new;
        at_old = at_new;

        % Compute the flow angle
        phi = atan(((1-an_old)*R)/(((1+at_old)*lambda*r)));   % [rad] : Flow angle

        % Compute local AoA alpha
        theta = theta_p + beta;         % [rad] : Twist angle = Global pitch angle + local twist    
        alpha = phi - theta;            % [rad] : Angle of attack = Flow angle - Twist angle

        % Read lift and drag coefficient from table look-up
        [Cl, Cd] = Airfoil_interpolation(n_airfoil, aoa_list, Cl_list, Cd_list, Thick_prof_list, k, alpha, itt);
        % disp([Cl,Cd, itt])
        


        % Compute Cn and Ct
        Cn = Cl*cos(phi) + Cd*sin(phi);     % [-] : Normal load coefficient
        Ct = Cl*sin(phi) - Cd*cos(phi);     % [-] : Tangential load coefficient

        % Compute CT
        sigma = (c*B)/(2*pi*r);                                 % [-] : Solidity factor
        F = 2/pi * acos(exp(-(B*(R-r))/(2*r*sin(abs(phi)))));   % [-] : Prandtl's tip loss corection
        CT = ((1-an_old)^2*Cn*sigma)/(sin(phi)^2);

        % Update an_new and at_new
        if an_old <= BetzLimit     % Momentum Theory Valid
            an_new = 1/((4*F*sin(phi)^2)/(sigma*Cn) + 1);
        else                    % Glauert Correction; Empirical
            beta_relaxation = 0.1;         % [-] : Relaxation factor needed (Given)
            a_star = (CT)/(4*F*(1-0.25*(5-3*an_old)*an_old));
            an_new = beta_relaxation*a_star+(1-beta_relaxation)*an_old; 
        end

        at_new = 1/((4*F*sin(phi)*cos(phi))/(sigma*Ct)-1);     % Assume no correction for a_t

        if (abs(an_new- an_old) > eps || abs(at_new-at_old) > eps)
            % Increase iteration step by one and keep on iterating
            % disp(abs(an_new- an_old))
            if alpha_mem < alpha
                alpha_mem = alpha;
            end
            n = n+1;
        else
            running = false;
        end
    end
    
    
    V_rel = V0*sqrt(((1+at_new)*lambda*r/R)^2 + (1-an_new)^2);   % [m/s] : Relative wind
    l = 0.5*rho*V_rel^2*c*Cl;                              % [N/m] : Lift per unit length
    d = 0.5*rho*V_rel^2*c*Cd;                              % [N/m] : Drag per unit length
    pn_r = l*cos(phi) + d*sin(phi);                        % [N/m] Normal load at radial position r
    pt_r = l*sin(phi) - d*cos(phi);                        % [N/m] Tangential load at radial position r
end
