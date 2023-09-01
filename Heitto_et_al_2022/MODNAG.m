%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for condensation model MODNAG                                                    %
% (Model for oligomerization and decomposition in nanoparticle growth)      %
%                                                                                                                   %
% By Arto Yli-Heitto, University of Eastern Finland, September 2020                %
% arto.heitto@uef.fi                                                                                        %
%                                                                                                                    %
% This routine calculates condensation of vapors on small monodisperse        %
% aerosol particles  and oligomerization and decomposition.                          %
%                                                                                                                    %
%                                                                                                                   %
% Vapor concentrations, RH and temperature are assumed to be constant     %
% during the simulation                                                                                 %
%                                                                                                                   %
% Modifications:                                                                                             %
% MABNAG-O build (T.Y. 2.7.2015)                                                                  %
%-------------------------------------------------------------------------             %
%                                                                                                                    %
% This is routine that is used when many model runs are made. This code       %
% makes the model calculations.                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version is for a system with 24+ condensing compounds.                  %
%                                                                                                                %
% The order of compounds for calculations is:                                              %
%   1 = water                                                                                               %
%   2 = ammonia                                                                                         %
%   3 = sulfuric acid                                                                                      %
%   4-10 = organics by VBS 2 - -4, no reactant                                              %
%   11-17 = organics by VBS 2 - -4, may oligomerize                                     %
%   18-24 = organics by VBS 2 - -4, may decompose                                     %
%   25-> = oligomerization and decomposing products                                %
%                                                                                                                  %
%                                                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear all
global R k NA whats watermark
global T press pvi RH
global M rhoi p_eq sigma D alpham
global mass_water mass_NH3 psat_w
global V_molec olig_ind decom_ind pro_ind P_sats
global activ product_numbers reversible
global simulation_compounds

input_file='input_MODNAG.xlsx'; %name of the excel file from where you read the input data

reversible=0; %if one, oligomerization reaction is reversible (i.e. decomposition happens to oligimerized compound)
watermark=1000; %whats the identifier for water in oligomerization and decomposition reactions

[data,txt,data_all]=xlsread(input_file); %read the data

org_ind=find(contains(txt,'organic'))-4; %Which are non reactant compounds?
olig_ind=find(contains(txt,'oligomerisative'))-4; %Which of the compounds can oligomerisate?
decom_ind=find(contains(txt,'decomposive'))-4; %Which of the compounds can decompose?

%-------------oligomerization reaction rate coefficient--------------------------
[RxnP_r,RxnP_c]=find(contains(txt,'Oligomerization'));
RxnP=data(RxnP_r-2:RxnP_r-3+sum(1:length(olig_ind)),RxnP_c-1:RxnP_c+3);
RxnP(RxnP(:,3)==0,:)=[];

%--------------------------------------------------------------------------

%-------------Decomposition reaction rate coefficient--------------------------
[DP_r,DP_c]=find(contains(txt,'Decomposition'));
if reversible ==0
DP=data(DP_r-2:DP_r-3+length(decom_ind),DP_c-1:DP_c+2);
else
    DP=data(DP_r-2:DP_r-3+length(decom_ind)+3,DP_c-1:DP_c+2);
end
DP(DP(:,2)==0,:)=[];

if reversible==0
    product_numbers=[RxnP(:,4);RxnP(:,5);DP(:,3);DP(:,4)];
    products_no_H2O=product_numbers;
    products_no_H2O(products_no_H2O==watermark)=[];
else
    product_numbers=[RxnP(:,4);RxnP(:,5)];
    products_no_H2O=product_numbers;
    products_no_H2O(products_no_H2O==watermark)=[];
end

simulation_compounds=[[1:org_ind(1)-1]'; org_ind; olig_ind; decom_ind; org_ind(1)-1+length(org_ind)+length(olig_ind)+length(decom_ind)+unique(products_no_H2O)];
%------------------Define basic parameters---------------------------------------
[BP_r,BP_c]=find(contains(txt,'BP'));

TotNumConc=data(BP_r-3,BP_c);
D_mean=data(BP_r-2,BP_c);
T=data(BP_r-1,BP_c); %Temperature (K)
sigma=data(BP_r,BP_c); %surface tension (N/m2)
RH=data(BP_r+1,BP_c); %Relative humidity (%)
press=data(BP_r+2,BP_c); %Pressure (Pa)
%--------------------------------------------------------------------------

%---------Physicochemical Parameters of components-------------------------
PCP=data(1:length(find(isfinite(data(:,5)))),1:length(find(isfinite(data(1,:))))+1);

rhoi=PCP(simulation_compounds,1);% Densities of the compounds (kg/m^3)
M=PCP(simulation_compounds,2)*1e-3;%molar mass [kg mol-1]
D=PCP(simulation_compounds,3);%gas diffusion coefficient[m2 s-1]
alpham=PCP(simulation_compounds,4);%surface accommodation coefficient on free substrate
C0=PCP(simulation_compounds,5);%effective saturation concentration[ug m-3]
dis1=PCP(simulation_compounds,6); %first dissociation constant for organic acids
dis2=PCP(simulation_compounds,7); %fsecond dissociation constant for organic diacids
whats=PCP(simulation_compounds,8); %What are the compounds
% -4 = HNO3
% -3 = H2SO4
% -2 = diacid
% -1 = monoacid
% 0 = water
% 3 = NH3
% 10 = neutral
Zg0=PCP(simulation_compounds,10);%initial gas phase concentration [cm-3]
VF0=PCP(simulation_compounds,11);%initial volume fraction of each component
pro_ind=find(contains(txt,'product'))-4;%Which of the compounds are products of oligomerization or decomposition?
pro_ind=pro_ind(unique(products_no_H2O));
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gas constant (J K^-1 mol^-1)
R = 8.314472;
%Boltzmann constant (m^2 kg s^-2 K^-1)
k = 1.3806503e-23;
%Avogadro constant (mol^-1)
NA = 6.0221415e+23;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: ambiet conditions and initial composition of a particle %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gas phase molecular concentrations #/m3
c_compounds = Zg0*1e+6;

% Saturation vapor pressure (Pa)
P_sats=C0.*1e-9*NA./M*k*T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: Properties of the compounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Converting concentrations (molecules/m^3) to pressures (Pa)
p_compounds = c_compounds.*k.*T;

%Calculating ambient partial pressure of water
%Coefficients from 1-component model file models.gro
A_w = [77.34491296e+0; 7235.424651e+0; 0.82e+1; 0.0057113e+0; 0.e+0];
psat_w = exp(A_w(1) - A_w(2)/T - A_w(3)*log(T) + A_w(4)*T + A_w(5)*T^2);
p_w = RH/100 * psat_w;

pvi = [p_w; p_compounds(2:end)];

%Initial radius of particle (m)
r0 = D_mean/2;

%Initial volume of particle (m^3)
V0 = 4/3*pi*r0^3; 

%Molecular volume (m^3)
V_molec = (M./NA)./rhoi;

%Molar volume (cm^3/mol)
V_molar = 1e6*M./rhoi;

%Initial composition of particles
Ni = VF0.*V0./V_molec;

% Masses of each component in a particle
mpi = Ni.*M ./NA;

% Mass of a particle
mp = sum(mpi);

% Density of the particle (kg/m^3)
rho = mp/V0;

%Initial radius of the particle
rp0 = (3*(mp/rho)/(4*pi)).^(1.0/3.0);

disp(['Initial radius of the particle is ' num2str(rp0) ' m.'])

% Kelvin effect
Ke = exp(2.*sigma.*V_molec.*NA./R./T./rp0);
%% 
% times
[time_r,time_c]=find(contains(txt,'Time grid'));  
time_data=data(time_r-3:time_r,time_c:length(find(isfinite(data(time_r,:)))));
starttimes=[];
for ii=1:size(time_data,2)
    if time_data(4,ii)==1 %number one here means, that times are logarhytmically spaced during this section
        starttimes=[starttimes logspace(time_data(1,ii), time_data(2,ii),time_data(3,ii))];
        if ii>1 && time_data(4,ii)==1
            starttimes(end-time_data(3,ii)+1)=[];
        end
    elseif time_data(4,ii)==2%number two here means, that times are linearly spaced during this section
        starttimes=[starttimes linspace(time_data(1,ii), time_data(2,ii),time_data(3,ii))];       
    end
end

 tout = [];
 output = [];
 jj = 0;
 time_start_loop = 0.;
 p_eq=zeros(length(whats),1);
 diss_frac=zeros(length(whats),1); %dissociated fraction of compounds

    % Inputs for the ode solver: only variables are the masses in the particulate phase
    input = mpi;

%%
    while time_start_loop < starttimes(end-2)
        jj = jj+1;
        start_time = log10(starttimes(jj));
        end_time = log10(starttimes(jj+1));
        time_start_loop = 10^start_time;
        time = logspace(start_time,end_time,10);

        % Masses of water and NH3 and equilibrium vapor pressures of organics
            [p_eq, input]=Calc_equil_water(input); 
      
            mass_water=input(whats==0);
            mass_NH3=input(whats==3);
            
     
        options = odeset('RelTol',1E-6,'AbsTol',1E-33);

        activ=p_eq./((input./M)./sum(input./M)); 
        activ(isinf(activ))=0;
        activ(isnan(activ))=0;
        [tout1, output1] = ode15s(@flx_eaim_MABNAGO_hd_pd_AH,time,input(3:end),options,whats,diss_frac,RxnP,DP);

        tout = [tout; tout1(2:end)];
        output = [output; repmat(mass_water,length(output1(2:end,1)),1) repmat(mass_NH3,length(output1(2:end,1)),1) output1(2:end,:)];
 
        input = output(end,:)'; 

        mass_all_1 = [mass_water mass_NH3 output1(1,:)];

        %Particle volume
        Vp = sum(mass_all_1'./M*NA.*V_molec);

        %Median radius:
        rp_median = median((3.*Vp./(4.*pi)).^(1/3));

        %If radius is > 40 nm the calculation will stop
        if rp_median > 40e-9
            break
        end
    end
%%
    clear mp Xmass Xmole rp

    mp = sum(output,2);
    Xmass = output./mp;

    n_times = size(output,1);

    Vp_out = sum(output./repmat(M',n_times,1)*NA.*repmat(V_molec',n_times,1),2);
    rp = (3.*Vp_out./(4.*pi)).^(1/3);

%--- Model output ---------------------------------------------------------
    model_result = [tout rp Xmass output];
%--------------------------------------------------------------------------
    
