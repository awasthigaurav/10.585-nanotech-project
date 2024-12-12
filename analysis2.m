length_vec = ["600", "800", "1000", "1200"];
colors = lines(4);
line_handles1 = gobjects(1, 4);
line_handles2 = gobjects(1, 4);

e = 1.6e-19;
Na = 6.022e23;
sigma = -0.06; % C/m2
kB = 1.38e-23;
T = 298.15;
eps = 78.54*8.8541e-12;
lB = e^2/(4*pi*eps*kB*T);
Ravg = 15e-9;
A = pi*Ravg^2;
Cratio = linspace(9,1100,15);
eta = 0.001; % viscosity
Cs = Cratio/2; % salt ions/m3
lambdaD = (8*pi*lB.*Cs).^(-0.5);
chi = 2*pi*lambdaD*lB*sigma/e;
alpha_ion = 0.8; % from emmerich et al.
b0 = 35e-9./(1+3e-4*(Cs/Na).^(2.4)); %converting Cs to mmol/l
betas = 50; % from emmerich et al.
Kosm = -0.5*2*sigma/(Ravg)*kB*T/(2*pi*eta*lB)*(1-asinh(chi)./chi+(1-alpha_ion)*b0./lambdaD.*(sqrt(1+chi.^2)-1)/(1+betas*abs(sigma)*lB^2/e));
    
for i = 1:4
    length = length_vec(i);
    L = str2double(length)*1e-9;
    data = readtable("Li_fig5a_" + length + "nm.csv", 'VariableNamingRule', 'preserve');
    x = data{:, 1}; 
    y = data{:, 2};
    data1 = readtable("Li_fig5d_" + length + "nm.csv", 'VariableNamingRule', 'preserve');
    x1 = data1{:, 1}; 
    y1 = data1{:, 2}; 
    
    Iosm = A/L.*Kosm.*log(Cratio);
    
    current_interp = interp1(log(Cratio), Iosm/1e-12, log(x), 'linear');
    SS_res = sum((y - current_interp).^2);        
    SS_tot = sum((y - mean(y)).^2);        
    R2 = 1 - (SS_res / SS_tot);             
    fprintf('R^2 for %.0f nm = %.4f\n', length_vec(i),R2);
    
    %% power
    mu = 4.8e11; %from siria et al, also similar by taking Li et al data
    G = 2*e^2*mu*Cs*A/L + e*mu*2*pi*Ravg/L*abs(sigma)*(1+1);
    Posm = Iosm.^2./(4*G);

    power_interp = interp1(log(Cratio), Posm/1e-12, log(x1), 'linear');
    SS_res = sum((y1 - power_interp).^2);        
    SS_tot = sum((y1 - mean(y1)).^2);        
    R2 = 1 - (SS_res / SS_tot);             
    fprintf('R^2 for P_{osm} = %.4f\n', R2);

    subplot(2,1,1)
    line_handles1(i) = semilogx(Cratio, Iosm/1e-12, 'Color', colors(i, :), 'LineWidth', 1.0);
    hold on
    semilogx(x, y, 'o', 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :))

    subplot(2,1,2)
    line_handles2(i) = semilogx(Cratio,Posm/1e-12, 'Color', colors(i, :), 'LineWidth', 1.0);
    hold on
    semilogx(x1, y1, 'o', 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :))
end

subplot(2,1, 1);
xlabel('C_{high}/C_{low}');
ylabel('I_{osm} (pA)');
title('Osmotic Current vs. Concentration Ratio');
legend(line_handles1,'600 nm', '800 nm', '1000 nm', '1200 nm', 'Location', 'best');

subplot(2,1, 2);
xlabel('C_{high}/C_{low}');
ylabel('P_{max} (pW)');
title('Max power vs. Concentration Ratio');
legend(line_handles2,'600 nm', '800 nm', '1000 nm', '1200 nm', 'Location', 'best');

%% replicating emmerich

e = 1.6e-19;
Na = 6.022e23;
sigma = -5; % C/m2
kB = 1.38e-23;
T = 298.15;
eps = 78.54*8.8541e-12;
lB = e^2/(4*pi*eps*kB*T);
w = 100e-9;
h = 3e-9;
A = w*h;
Cratio = linspace(3,1100,15);
eta = 0.001; % viscosity
Cs = Cratio/2; % salt ions/m3
lambdaD = (8*pi*lB.*Cs).^(-0.5);
chi = 2*pi*lambdaD*lB*sigma/e;
alpha_ion = 0.8; % from emmerich et al.
b0 = 35e-9./(1+3e-4*(Cs/Na).^(2.4)); %converting Cs to mmol/l
betas = 50; % from emmerich et al.
Kosm = -2*sigma/(h)*kB*T/(2*pi*eta*lB)*(1-asinh(chi)./chi+(1-alpha_ion)*b0./lambdaD.*(sqrt(1+chi.^2)-1)/(1+betas*abs(sigma)*lB^2/e));
L = 10000e-9;
Iosm = A/L.*Kosm.*log(Cratio);

data = readtable("emmerich_fig4a_replication.xlsx",'VariableNamingRule', 'preserve');
x = data{:, 1}; 
y = data{:, 2}*1e3;
current_interp = interp1(log(Cratio), Iosm/1e-12, log(x), 'linear');
SS_res = sum((y - current_interp).^2);        
SS_tot = sum((y - mean(y)).^2);        
R2 = 1 - (SS_res / SS_tot);             
fprintf('R^2 = %.4f\n',R2);

figure()
semilogx(Cratio,Iosm/1e-12,'LineWidth',1.5)
hold on
semilogx(x, y, 'o','MarkerEdgeColor','black','MarkerFaceColor','black')
hold off
legend('Model prediction','Data','Location','best')
xlabel('C_{high}/C_{low}');
ylabel('I_{osm} (pA)');
title('Osmotic Current vs. Concentration Ratio');