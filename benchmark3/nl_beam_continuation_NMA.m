%========================================================================
% DESCRIPTION: 
% Investigation of a curved beam considering geometric nonlinearities.
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.1 Copyright (C) 2019  Malte Krack  
%										(malte.krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(johann.gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
clear; clc;
close all;
addpath('00_SRC');
addpath('00_SRC/MechanicalSystems');
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 

%% Define system

% Fundamental parameters
Dmod = [.38 .09 .08]*.01; % first, third, fifth bending modes of flat beam
Nmod = 3;
setup = 'New_Design_Steel';

thickness = .001;
R=3;

[L,rho,E,~,PHI,~,gam] = beams_for_everyone(setup,Nmod*2-1,thickness);
PHI_L_2 = PHI(L/2);
PHI_L_2 = PHI_L_2(1:2:end);
gam = gam(1:2:end);

% load nonlinear coefficients (can be found e.g. via IC-method)
load(['beam_msh_80_4_1_3t_steel_' num2str(thickness*1000) 'mm_R' num2str(R) 'm.mat'])
om = model.omega;

% Properties of the underlying linear system
M = eye(Nmod);
D = diag(2*Dmod(:).*om(:));
K = diag(om.^2);

% polynomial terms
p_quad = zeros(sum(1:Nmod),Nmod);
p_cub = zeros(sum(cumsum(1:Nmod)),Nmod);
ctr_quad = 1; ctr_cub = 1;

for jj = 1:Nmod
    for kk = jj:Nmod
        % quadratic terms
        p_quad(ctr_quad,jj) = p_quad(ctr_quad,jj)+1;
        p_quad(ctr_quad,kk) = p_quad(ctr_quad,kk)+1;
        ctr_quad = ctr_quad+1;
        for ll = kk:Nmod
            % cubic terms
            p_cub(ctr_cub,jj) = p_cub(ctr_cub,jj)+1;
            p_cub(ctr_cub,kk) = p_cub(ctr_cub,kk)+1;
            p_cub(ctr_cub,ll) = p_cub(ctr_cub,ll)+1;
            ctr_cub = ctr_cub+1;
        end
    end
end

p = [p_quad; p_cub];

% coefficients
E=zeros(sum(cumsum(1:Nmod)),Nmod);

for rr = 1:Nmod
    ctr = 1;
        for jj = 1:Nmod
            for kk = jj:Nmod
                % quadratic coeffs
                E(ctr,rr) = model.a(jj,kk,rr);
                ctr = ctr+1;
            end
        end
%         ctr = 1;
    for jj = 1:Nmod
        for kk = jj:Nmod
            for ll = kk:Nmod
                % cubic coeffs
                E(ctr,rr) = model.b(jj,kk,ll,rr);
                ctr = ctr+1;
            end
        end
    end
end

% Fundamental harmonic of external forcing
Fex1 = gam;

% Define oscillator as system with polynomial stiffness nonlinearities
oscillator = System_with_PolynomialStiffnessNonlinearity(M,D,K,p,E,Fex1);

% Number of degrees of freedom
n = oscillator.n;

%% Compute frequency response using harmonic balance
analysis = 'FRF';
H = 20;              % harmonic order
N=2*3*H+1;
Ntd = 1e3;

% Analysis parameters
Om_e = 0;      % start frequency
Om_s = 3*om(1);     % end frequency

% Excitation levels
exc_lev = [30];
Om = cell(size(exc_lev));
for iex=1:length(exc_lev)
    % Set excitation level
    oscillator.Fex1 = Fex1*exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*M + 1i*Om_s*D + K)\oscillator.Fex1;
    y0 = zeros((2*H+1)*length(Q1),1);
    y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];
    qscl = max(abs((-om(1)^2*M + 1i*om(1)*D + K)\oscillator.Fex1));
    
    % Solve and continue w.r.t. Om
    ds = 50; % -> better for exc_lev = 50
        
    Dscale = [1e-6*ones(length(y0),1);(Om_s+Om_e)/2];
    Sopt = struct('Dscale',Dscale,'dynamicDscale',1,'jac','full','stepmax',1e4);
    X = solve_and_continue(y0,...
        @(X) HB_residual(X,oscillator,H,N,analysis),...
        Om_s,Om_e,ds,Sopt);
    
    % Interpret solver output
    Om{iex} = X(end,:);
    Q_HB = X(1:end-1,:);
           
    % Define amplitude as magnitude of the fundamental harmonic of the
    % first coordinate's displacement
    tau = linspace(0,2*pi,Ntd+1); tau = tau(1:end-1);    
    w_L_2_sum = zeros(Ntd,length(Om{1}));                                   % lines = time steps, columns = continuation points
    for k = 1:n
        Qc{k} = [Q_HB(k,:);Q_HB(n+k:2*n:end,:)-1i*Q_HB(2*n+k:2*n:end,:)];   % get complex coefficients
        w_L_2{k} = PHI_L_2(k)*real(exp(1i*tau(:)*(0:H))*Qc{k});             % get displacement at center caused by each mode in time domain
        w_L_2_sum = [w_L_2_sum + w_L_2{k}];                                 % sum up to actual displacement at center in time domain
        
        q{k} = real(exp(1i*tau(:)*(0:H))*Qc{k});
        a_q(k,:) = (max((q{k}))-min((q{k})))/2;
    end
    
    a_w_L_2= (max((w_L_2_sum))-min((w_L_2_sum)))/2;                         % compute peak to peak amplitude
    
end

% Illustrate frequency response
figure; hold on; box on
plot(Om{1}/2/pi,a_w_L_2*1000,'k-','linewidth',1.5);
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('$\hat{w}_{L/2}$ in mm');
title(['thickness = ' num2str(thickness*1000) 'mm, $R=$' num2str(R) 'm, $\hat{\ddot{a}}_0=$' num2str(exc_lev) 'm/s$^2$']);
    
%% NMA
H=7;
N=2*3*H+1;

% Linear modal analysis
[PHI_lin,OM2] = eig(oscillator.K,oscillator.M);
[om_lin,ind] = sort(sqrt(diag(OM2)));
PHI_lin = PHI_lin(:,ind);

analysis = 'NMA';

imod = 1;           % mode to be analyzed
log10a_s = -7;    % start vibration level (log10 of modal mass)
log10a_e = -3.2;       % end vibration level (log10 of modal mass)
inorm = 1;          % coordinate for phase normalization

% Initial guess vector y0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_lin(imod); phi = PHI_lin(:,imod);
Psi = zeros((2*H+1)*Nmod,1);
Psi(Nmod+(1:Nmod)) = phi;
x0 = [Psi;om;0];

ds      = .1;
Sopt    = struct('Dscale',[1e-6*ones(size(x0,1)-2,1);1;1e-1;1],...
    'dynamicDscale',1,'stepmax',5e4);
% Sopt    = struct('stepmax',5e4);
[X_HB,Solinfo,Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H,N,analysis,inorm),...
    log10a_s,log10a_e,ds, Sopt);

Psi_HB = X_HB(1:end-3,:);
om_HB = X_HB(end-2,:);
del_HB = X_HB(end-1,:);
log10a_NMA = X_HB(end,:);
a_NMA = 10.^log10a_NMA;
Q_HB = Psi_HB.*repmat(a_NMA,size(Psi_HB,1),1);

% Define amplitude as magnitude of the fundamental harmonic of the
% first coordinate's displacement
tau = linspace(0,2*pi,Ntd+1); tau = tau(1:end-1);
w_L_2_NMA_sum = zeros(Ntd,length(om_HB));                                   % lines = time steps, columns = continuation points
for k = 1:n
    Qc_NMA{k} = [Q_HB(k,:);Q_HB(n+k:2*n:end,:)-1i*Q_HB(2*n+k:2*n:end,:)];   % get complex coefficients
    w_L_2_NMA{k} = PHI_L_2(k)*real(exp(1i*tau(:)*(0:H))*Qc_NMA{k});         % get displacement at center caused by each mode in time domain
    w_L_2_NMA_sum = [w_L_2_NMA_sum + w_L_2_NMA{k}];                         % sum up to actual displacement at center in time domain
    
    q_NMA{k} = real(exp(1i*tau(:)*(0:H))*Qc_NMA{k});
    a_q_NMA(k,:) = (max((q_NMA{k}))-min((q_NMA{k})))/2;
end

a_w_L_2_NMA= (max((w_L_2_NMA_sum))-min((w_L_2_NMA_sum)))/2;                 % compute peak to peak amplitude

figure(1)
hold on, box on
plot(om_HB/2/pi,1000*abs(a_w_L_2_NMA),'color',[1 .2 .3],'linewidth',1.5)
plot(Om{1}/2/pi,a_w_L_2*1000,'k-','linewidth',1.5);
xlim([260 420])
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('$\hat{w}_{L/2}$ in mm');
title(['thickness = ' num2str(thickness*1000) 'mm, $R=$' num2str(R) 'm, $\hat{\ddot{a}}_0=$' num2str(exc_lev) 'm/s$^2$']);

figure;
plot(om_HB/2/pi,a_q_NMA,'linewidth',1.5)
legend('$q_1$ (first bending)','$q_2$ (third bending)','$q_3$ (fifth bending)')
xlabel('$\omega$ in Hz');
ylabel('modal coordinates $q$ of NMA')
xlim([min(om_HB/2/pi) max(om_HB/2/pi)]);
title(['thickness = ' num2str(thickness*1000) 'mm, $R=$' num2str(R) 'm, $\hat{\ddot{a}}_0=$' num2str(exc_lev) 'm/s$^2$']);