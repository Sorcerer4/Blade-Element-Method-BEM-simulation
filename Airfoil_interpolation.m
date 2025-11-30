%% Airfoil Properties Interpolation 
% To be done for each thichness k and angle of attack alpha

function [Cl, Cd] = Airfoil_interpolation(n_airfoil, aoa_list, Cl_list, Cd_list, Thick_prof_list, k, alpha, itt)
    
    Cl_thick = zeros(n_airfoil,1);
    Cd_thick = zeros(n_airfoil,1);

    disp([rad2deg(alpha), itt])

    % Interpolate the values to the different thicknesses
    for s=1:n_airfoil       % s indicate the airfoil section
        Cl_thick(s,1)=interp1(aoa_list(:,s),Cl_list(:,s),alpha);
        Cd_thick(s,1)=interp1(aoa_list(:,s),Cd_list(:,s),alpha);
    end
    
    % then interpolate to the actual thickness
    Cl = interp1(Thick_prof_list(:),Cl_thick(:),k);
    Cd = interp1(Thick_prof_list(:),Cd_thick(:),k);

end