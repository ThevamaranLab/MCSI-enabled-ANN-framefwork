                                  clc
clear all
close all
set(0,'defaultfigurewindowstyle','docked')
warning off
set(groot,'defaultAxesFontSize',18)
set(0,'defaultfigurecolor',[1 1 1])
set(0, 'DefaultLineLineWidth', 3);

%% Non-Architected

load('Non_Architected_stress_strain.mat')

for ind=1:length(Na_Strain_L)

    xv = Na_Strain_L{ind}(100:end);
    yv = Na_Stress_L{ind}(100:end);
    f = fit(xv,yv,'El*x^ll');
    Na_samp_dat(ind,10) = f.El;
    Na_samp_dat(ind,11) = f.ll;

    xv2 = Na_Strain_UL{ind}(1:end-100);
    yv2 = Na_Stress_UL{ind}(1:end-100);
    f1 = fit(xv2,yv2,'Eul*x^lul');
    Na_samp_dat(ind,12) = f1.Eul;
    Na_samp_dat(ind,13) = f1.lul;

    Na_samp_dat(ind,14) = (f1.lul-f.ll)/(f1.lul+1);  % Damping capacity from fit
    Na_samp_dat(ind,15) = Na_samp_dat(ind,7)/2.26; % Relative density
    
    Eps_c = 1-1.733*Na_samp_dat(ind,15); % strain at the onset to densification
    Na_samp_dat(ind,16) = (f.ll+1)./(Eps_c);  % hbar_cr
    E_l_bar(ind) = Na_samp_dat(ind,10)/15000;  % Relative modulus

    Na_samp_dat(ind,17) = 1/(((Eps_c)^(Na_samp_dat(ind,11)))*E_l_bar(ind)); % Abar_cr
    
    Na_samp_dat(ind,18) = trapz(Na_Strain_L{ind},Na_Stress_L{ind}); % Loading energy in MJ/m^3
    Na_samp_dat(ind,19) = -trapz(Na_Strain_UL{ind},Na_Stress_UL{ind}); % unloading area MJ/m^3
    Na_samp_dat(ind,20) = Na_samp_dat(ind,18)-Na_samp_dat(ind,19); % Dissipated area MJ/m^3
    Na_samp_dat(ind,21) = 0.5*(Na_samp_dat(ind,18)+Na_samp_dat(ind,19))*2*(1/(0.5)^2); % average modulus in MPa
    Na_samp_dat(ind,22) = Na_samp_dat(ind,20)/Na_samp_dat(ind,15); % Specific Dissipated area MJ/m^3
    Na_samp_dat(ind,23) = Na_samp_dat(ind,21)/Na_samp_dat(ind,15); % Specific average modulus in MPa
    

    figure
    plot(xv,yv,'LineWidth',3,'Color','k')
    hold on
    plot(xv,Na_samp_dat(ind,10)*(xv.^Na_samp_dat(ind,11)),'LineWidth',3,'Color','k','LineStyle','--')   
    plot(xv2,yv2,'LineWidth',3,'Color','r')
    plot(xv,Na_samp_dat(ind,12)*(xv.^Na_samp_dat(ind,13)),'LineWidth',3,'Color','r','LineStyle','--')

end
%}
%% Nested

load('Nested_stress_strain.mat')

for ind=1:length(Ns_Strain_L)

    xv = Ns_Strain_L{ind}(100:end);
    yv = Ns_Stress_L{ind}(100:end);
    f = fit(xv,yv,'El*x^ll');
    Ns_samp_dat(ind,10) = f.El;
    Ns_samp_dat(ind,11) = f.ll;

    xv2 = Ns_Strain_UL{ind}(1:end-100);
    yv2 = Ns_Stress_UL{ind}(1:end-100);
    f1 = fit(xv2,yv2,'Eul*x^lul');
    Ns_samp_dat(ind,12) = f1.Eul;
    Ns_samp_dat(ind,13) = f1.lul;

    Ns_samp_dat(ind,14) = (f1.lul-f.ll)/(f1.lul+1);  % Damping capacity from fit
    Ns_samp_dat(ind,15) = Ns_samp_dat(ind,7)/2.26; % Relative density
    
    Eps_c = 1-2.287*Ns_samp_dat(ind,15); % strain at the onset to densification
    Ns_samp_dat(ind,16) = (f.ll+1)./(Eps_c);  % hbar_cr
    E_l_bar(ind) = Ns_samp_dat(ind,10)/15000;  % Relative modulus

    Ns_samp_dat(ind,17) = 1/(((Eps_c)^(Ns_samp_dat(ind,11)))*E_l_bar(ind)); % Abar_cr

    Ns_samp_dat(ind,18) = trapz(Ns_Strain_L{ind},Ns_Stress_L{ind}); % Loading energy in MJ/m^3
    Ns_samp_dat(ind,19) = -trapz(Ns_Strain_UL{ind},Ns_Stress_UL{ind}); % unloading area MJ/m^3
    Ns_samp_dat(ind,20) = Ns_samp_dat(ind,18)-Ns_samp_dat(ind,19); % Dissipated area MJ/m^3
    Ns_samp_dat(ind,21) = 0.5*(Ns_samp_dat(ind,18)+Ns_samp_dat(ind,19))*2*(1/(0.5)^2); % average modulus in MPa
    Ns_samp_dat(ind,22) = Ns_samp_dat(ind,20)/Ns_samp_dat(ind,15); % Specific Dissipated area MJ/m^3
    Ns_samp_dat(ind,23) = Ns_samp_dat(ind,21)/Ns_samp_dat(ind,15); % Specific average modulus in MPa
    

    figure
    plot(xv,yv,'LineWidth',3,'Color','k')
    hold on
    plot(xv,Ns_samp_dat(ind,10)*(xv.^Ns_samp_dat(ind,11)),'LineWidth',3,'Color','k','LineStyle','--')   
    plot(xv2,yv2,'LineWidth',3,'Color','r')
    plot(xv,Ns_samp_dat(ind,12)*(xv.^Ns_samp_dat(ind,13)),'LineWidth',3,'Color','r','LineStyle','--')

end

%% Cylindrical

load('cylindrical_stress_strain.mat')

for ind=1:length(Cy_Strain_L)

    xv = Cy_Strain_L{ind}(100:end);
    yv = Cy_Stress_L{ind}(100:end);
    f = fit(xv,yv,'El*x^ll');
    Cy_samp_dat(ind,10) = f.El;
    Cy_samp_dat(ind,11) = f.ll;

    xv2 = Cy_Strain_UL{ind}(1:end-100);
    yv2 = Cy_Stress_UL{ind}(1:end-100);
    n_i = find(yv2 < 0);

    if isempty(n_i)==1
    f1 = fit(xv2,yv2,'Eul*x^lul');
    figure
    plot(xv2,yv2,'LineWidth',3,'Color','r')
    else
    xv3 = xv2(1:n_i(1)-1);
    yv3 = yv2(1:n_i(1)-1);
    f1 = fit(xv3,yv3,'Eul*x^lul');
    figure
    plot(xv3,yv3,'LineWidth',3,'Color','r')
    end
    Cy_samp_dat(ind,12) = f1.Eul;
    Cy_samp_dat(ind,13) = f1.lul;

    Cy_samp_dat(ind,14) = (f1.lul-f.ll)/(f1.lul+1);  % Damping capacity from fit
    Cy_samp_dat(ind,15) = Cy_samp_dat(ind,7)/2.26; % Relative density
    
    Eps_c = 1-2.287*Cy_samp_dat(ind,15); % strain at the onset to densification
    Cy_samp_dat(ind,16) = (f.ll+1)./(Eps_c);  % hbar_cr
    E_l_bar(ind) = Cy_samp_dat(ind,10)/15000;  % Relative modulus

    Cy_samp_dat(ind,17) = 1/(((Eps_c)^(Cy_samp_dat(ind,11)))*E_l_bar(ind)); % Abar_cr

    Cy_samp_dat(ind,18) = trapz(Cy_Strain_L{ind},Cy_Stress_L{ind}); % Loading energy in MJ/m^3
    Cy_samp_dat(ind,19) = -trapz(Cy_Strain_UL{ind},Cy_Stress_UL{ind}); % unloading area MJ/m^3
    Cy_samp_dat(ind,20) = Cy_samp_dat(ind,18)-Cy_samp_dat(ind,19); % Dissipated area MJ/m^3
    Cy_samp_dat(ind,21) = 0.5*(Cy_samp_dat(ind,18)+Cy_samp_dat(ind,19))*2*(1/(0.5)^2); % average modulus in MPa
    Cy_samp_dat(ind,22) = Cy_samp_dat(ind,20)/Cy_samp_dat(ind,15); % Specific Dissipated area MJ/m^3
    Cy_samp_dat(ind,23) = Cy_samp_dat(ind,21)/Cy_samp_dat(ind,15); % Specific average modulus in MPa
    

    hold on
    plot(xv,Cy_samp_dat(ind,12)*(xv.^Cy_samp_dat(ind,13)),'LineWidth',3,'Color','r','LineStyle','--')
    plot(xv,yv,'LineWidth',3,'Color','k')
    plot(xv,Cy_samp_dat(ind,10)*(xv.^Cy_samp_dat(ind,11)),'LineWidth',3,'Color','k','LineStyle','--')   

end

%% Concentric

load('Concentric_stress_strain.mat')

for ind=1:length(Co_Strain_L)

    xv = Co_Strain_L{ind}(100:end);
    yv = Co_Stress_L{ind}(100:end);
    f = fit(xv,yv,'El*x^ll');
    Co_samp_dat(ind,10) = f.El;
    Co_samp_dat(ind,11) = f.ll;

    xv2 = Co_Strain_UL{ind}(1:end-100);
    yv2 = Co_Stress_UL{ind}(1:end-100);
    n_i = find(yv2 < 0);

    if isempty(n_i)==1
    f1 = fit(xv2,yv2,'Eul*x^lul');
    figure
    plot(xv2,yv2,'LineWidth',3,'Color','r')
    else
    xv3 = xv2(1:n_i(1)-1);
    yv3 = yv2(1:n_i(1)-1);
    f1 = fit(xv3,yv3,'Eul*x^lul');  
    figure
    plot(xv3,yv3,'LineWidth',3,'Color','r')
    end
    Co_samp_dat(ind,12) = f1.Eul;
    Co_samp_dat(ind,13) = f1.lul;

    Co_samp_dat(ind,14) = (f1.lul-f.ll)/(f1.lul+1);  % Damping capacity from fit
    Co_samp_dat(ind,15) = Co_samp_dat(ind,7)/2.26; % Relative density
    
    Eps_c = 1-2.287*Co_samp_dat(ind,15); % strain at the onset to densification
    Co_samp_dat(ind,16) = (f.ll+1)./(Eps_c);  % hbar_cr
    E_l_bar(ind) = Co_samp_dat(ind,10)/15000;  % Relative modulus

    Co_samp_dat(ind,17) = 1/(((Eps_c)^(Co_samp_dat(ind,11)))*E_l_bar(ind)); % Abar_cr

    Co_samp_dat(ind,18) = trapz(Co_Strain_L{ind},Co_Stress_L{ind}); % Loading energy in MJ/m^3
    Co_samp_dat(ind,19) = -trapz(Co_Strain_UL{ind},Co_Stress_UL{ind}); % unloading area MJ/m^3
    Co_samp_dat(ind,20) = Co_samp_dat(ind,18)-Co_samp_dat(ind,19); % Dissipated area MJ/m^3
    Co_samp_dat(ind,21) = 0.5*(Co_samp_dat(ind,18)+Co_samp_dat(ind,19))*2*(1/(0.5)^2); % average modulus in MPa
    Co_samp_dat(ind,22) = Co_samp_dat(ind,20)/Co_samp_dat(ind,15); % Specific Dissipated area MJ/m^3
    Co_samp_dat(ind,23) = Co_samp_dat(ind,21)/Co_samp_dat(ind,15); % Specific average modulus in MPa
    

    hold on
    plot(xv,Co_samp_dat(ind,12)*(xv.^Co_samp_dat(ind,13)),'LineWidth',3,'Color','r','LineStyle','--')
    plot(xv,yv,'LineWidth',3,'Color','k')
    plot(xv,Co_samp_dat(ind,10)*(xv.^Co_samp_dat(ind,11)),'LineWidth',3,'Color','k','LineStyle','--')    

end

%%


