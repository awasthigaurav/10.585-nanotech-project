Na = 6.022e23;
e = 1.6e-19;
L = 9600e-9*100; % ultra-long CNT, approx 1 mm
D = 2e-9;
V = 1e-3*(150:50:1000);
E = V/L;
muH = 3.62e-7*(1+39*exp(-(D*1e9-0.81)/0.33));
sigma = 0.5; % for activated channels, taken from emmerich
CHclosed = sigma*2*(pi*D)/(e);
CHopen = CHclosed*1.2;

eps = 500; % suppl fig s12 choi et al
eta = 0.01;
nub = sqrt(e*CHclosed*muH*E.^2*eps/(eta*8*pi));
nu0 = e*CHopen*(muH*E+nub).*E*eps./(eta*8*pi*nub);
Iopen = e*E*muH*CHopen+e*nu0*CHopen;
Iclosed = e*nub*CHclosed;

figure(1)
plot(V,(Iopen-Iclosed)/1e-12,'LineWidth',1)
hold on
data = readtable("choi_replication.csv",'VariableNamingRule', 'preserve');
data_voltage = 1e-3*data{:, 1}; %converting mV to V
data_current = data{:,2}; % in pA
plot(data_voltage,data_current,'o','MarkerEdgeColor','black','MarkerFaceColor','black')
hold off
legend('Model prediction','Data','Location','best')
xlabel('Voltage (V)')
ylabel('Current (pA)')
xlim([0.19 1.05])
title('Current vs Voltage')
%%  emmerich 2e - activated
w = 100e-9;
h = 3e-9;
L = 10000e-9;
V = -0.1:0.01:0.1;
E = V/L;
D = 4*w*h/(2*(w+h));
muH = 3.62e-7*(1+39*exp(-(D*1e9-0.81)/0.33));
eps = 100; %suppl fig s12 choi et al
eta = 0.01;
conc_vec = [1,10,100,300];
sigma_vec = [0.3,0.5,1,2];
colors = lines(4);
line_handles1 = gobjects(1, 4);

for i = 1:4
    data = readtable("emmerich_fig2e.xlsx",'VariableNamingRule', 'preserve');
    emmerich_current =  data{:, i}'; 
    sigma = sigma_vec(i); % estimated from emmerich suppl table 4.1, this is where the conc dependence comes from
    CHclosed = sigma*2*(pi*D)/(e);
    CHopen = CHclosed*1.2;

    nub = sqrt(e*CHclosed*muH*E.^2*eps/(eta*8*pi)).*sign(V);
    nu0 = e*CHopen*(muH*E+nub).*E*eps./(eta*8*pi*nub);
    Iopen = e*E*muH*CHopen+e*nu0*CHopen;
    Iclosed = e*nub*CHclosed;
    
    % figure(2)
    line_handles1(i) = plot(V,Iclosed/1e-9,'LineWidth',1.5,'Color', colors(i, :));
    hold on
    plot(V,emmerich_current,'o','MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :))
    hold on

    SS_res = sum((Iclosed/1e-9 - emmerich_current).^2);        
    SS_tot = sum((emmerich_current - mean(emmerich_current)).^2);        
    R2 = 1 - (SS_res / SS_tot);
    fprintf('R^2 for %.0f mM = %.4f\n', conc_vec(i),R2);
end
hold off
legend(line_handles1,'1mM','10mM','100mM','300mM','Location','best')
xlabel('Voltage (V)')
ylabel('Current (nA)')
title('Current vs voltage - activated channel')

%% emmerich 2a - pristine
w = 120e-9;
h = 3e-9;
L = 9600e-9;
V = -0.1:0.01:0.1;
E = V/L;
D = 4*w*h/(2*(w+h));
muH = 3.62e-7*(1+39*exp(-(D*1e9-0.81)/0.33));
eps = 100; %suppl fig s12 choi et al
eta = 0.01;
conc_vec = [1,10,100,300];
sigma_vec = [2e-3,1e-2,5.4e-2,0.15];
colors = lines(length(conc_vec));
line_handles1 = gobjects(1, length(conc_vec));

for i = 1:length(conc_vec)
    data = readtable("emmerich_fig2a.xlsx",'VariableNamingRule', 'preserve');
    emmerich_current =  data{:, 2*i}';
    emmerich_current = emmerich_current(~isnan(emmerich_current));
    emmerich_voltage = data{:,2*i-1}';
    emmerich_voltage = emmerich_voltage(~isnan(emmerich_voltage));
    sigma = sigma_vec(i); % estimated from emmerich suppl table 4.1, this is where the conc dependence comes from
    CHclosed = sigma*2*(pi*D)/(e);
    CHopen = CHclosed*1.2;

    nub = sqrt(e*CHclosed*muH*E.^2*eps/(eta*8*pi)).*sign(V);
    nu0 = e*CHopen*(muH*E+nub).*E*eps./(eta*8*pi*nub);
    Iopen = e*E*muH*CHopen+e*nu0*CHopen;
    Iclosed = e*nub*CHclosed;
    
    % figure(2)
    line_handles1(i) = plot(V,Iclosed/1e-9,'LineWidth',1.5,'Color', colors(i, :));
    hold on
    plot(emmerich_voltage,emmerich_current,'o','MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :))
    hold on

    current_interp = interp1(V, Iclosed/1e-9, emmerich_voltage, 'linear');
    SS_res = sum((emmerich_current - current_interp).^2);        
    SS_tot = sum((emmerich_current - mean(emmerich_current)).^2);        
    R2 = 1 - (SS_res / SS_tot);                      

    fprintf('R^2 for %.0f mM = %.4f\n', conc_vec(i),R2);
end
hold off
legend(line_handles1,'1mM','10mM','100mM','300mM','Location','best')
xlabel('Voltage (V)')
ylabel('Current (nA)')
title('Current vs voltage - pristine channel')

%% changing length, cross-sectional area and voltage to match Choi et al. while still maintaining rectangular cross-section
Na = 6.022e23;
e = 1.6e-19;
w = 100e-9;
h = 1e-9;
L = 9600e-9*100;
V = 1e-3*(200:50:1000);%-0.1:0.01:0.1;
E = V/L;
D = 4*w*h/(2*(w+h));
muH = 3.62e-7*(1+39*exp(-(D*1e9-0.81)/0.33));
sigma = 0.5; % for actiavted channels
CHclosed = sigma*2*(pi*D)/(e);
CHopen = CHclosed*1.2;

eps = 500;
eta = 0.01;
nub = sqrt(e*CHclosed*muH*E.^2*eps/(eta*8*pi)).*sign(V);
nu0 = e*CHopen*(muH*E+nub).*E*eps./(eta*8*pi*nub);
Iopen = e*E*muH*CHopen+e*nu0*CHopen;
Iclosed = e*nub*CHclosed;

figure(2)
plot(V,(Iopen-Iclosed)/1e-12,'LineWidth',1.5)
hold on


data = readtable("choi_replication.csv",'VariableNamingRule', 'preserve');
data_voltage = 1e-3*data{:, 1}; %converting mV to V
data_current = data{:,2}; % in pA
plot(data_voltage,data_current,'o','MarkerEdgeColor','black','MarkerFaceColor','black')
hold off
legend('Model prediction','Data','Location','best')
xlabel('Voltage (V)')
ylabel('Current (pA)')
title('Hypothetical ultra-long rectangular nanochannel')
%% accounting for entrance effects
Na = 6.022e23;
e = 1.6e-19;
w = 100e-9;
h = 3e-9;
L = 10000e-9;
V = -0.1:0.01:0.1;
E = V/L;
D = 4*w*h/(2*(w+h));
muH = 3.62e-7*(1+39*exp(-(D*1e9-0.81)/0.33));
eps = 100; %suppl fig s12 choi et al
eta = 0.01;
conc_vec = [1,10,100,300];
sigma_vec = [0.3,0.5,1,2];
colors = lines(4);
line_handles1 = gobjects(1, 4);

for i = 1:4
    data = readtable("emmerich_fig2e.xlsx",'VariableNamingRule', 'preserve');
    emmerich_current =  data{:, i}'; 
    sigma = sigma_vec(i); % estimated from emmerich suppl table 4.1, this is where the conc dependence comes from
    CHclosed = sigma*2*(pi*D)/(e);
    CHopen = CHclosed*1.2;

    nub = sqrt(e*CHclosed*muH*E.^2*eps/(eta*8*pi)).*sign(V);
    nu0 = e*CHopen*(muH*E+nub).*E*eps./(eta*8*pi*nub);
    Iopen = e*E*muH*CHopen+e*nu0*CHopen;
    Iclosed = e*nub*CHclosed;
    
    figure()
    line_handles1(i) = plot(V,Iclosed/1e-9,'LineWidth',1.0,'Color', colors(i, :));
    hold on
    plot(V,emmerich_current,'o','MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :))
    hold on
    Raccess = 0.2/(4*D/2)*500/conc_vec(i); % seawater has around 0.5 M NaCl and a resistivity of 0.2 ohm-m, rescaling with that
    Iaccess = V./(V./Iclosed+Raccess);
    plot(V,Iaccess/1e-9,'LineStyle','--','LineWidth',1.5,'Color', colors(i, :));
    hold off
    legend('Model','Data','Model with R_{access}','Location','best')
    title(['Accounting for R_{access}; C=',num2str(conc_vec(i)),'mM'])
    xlabel('Voltage (V)')
    ylabel('Current (nA)')
end
hold off

