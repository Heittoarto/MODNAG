function flx = flx_MODNAG(time, input,whats,D_frac,RxnP,DP) 

global R k NA
global T press pvi olig_ind decom_ind pro_ind
global M rhoi D alpham 
global mass_water mass_NH3
global V_molec product_numbers
global activ watermark reversible
global simulation_compounds

%Masses of each component
mpi=zeros(length(whats),1);
mpi(whats==0)=mass_water;
mpi(whats==3)=mass_NH3;
mpi(whats<0 | whats>=10)=input;

%Particle mass
mp = sum(mpi);

%Mass fractions
Xmass = mpi./mp;

%Mole fractions
Xmole = (Xmass./M)./sum(Xmass./M);

%Particle radius
Vp = sum(mpi./M.*NA.*V_molec);
rp = (3.0*Vp/(4.0*pi)).^(1.0/3.0);

dp = 2.0*rp;

%--- Transition regime corrections ----------------------------------------
%According to Lehtinen & Kulmala and Neiminen et al.

%Calculating diffusion coefficient of particle
Ts=120; %Sutherland's coefficient for water (Licht and Stechert, 1944)
T0=291.15; %Reference temperaure (K)
visc_air0=18.27e-6; %Gas viscosity of air in reference temperature (kg/m/s)
visc_air=visc_air0*(T0+Ts)/(T+Ts)*(T/T0)^(3/2); %Gas viscosity of air (kg/m/s)
M_air = 29.0e-3;  %Molar mass of air (kg/mol)
lambda_air = 2*visc_air/(press*(8*M_air/(pi*R*T))^(1/2));  %Mean free bath of air molecules
Kn_air = 2*lambda_air/dp;  %Knudsen number using air mean free bath
Cc = 1+Kn_air*(1.257+0.4*exp(-1.1/Kn_air)); %Cunningham slip correction factor
Diffp = k*T*Cc/(3*pi*visc_air*dp);  %Diffusion coefficient of particle

%Mass of vapor molecules
mv = M./NA;

%Mean thermal speed
cv = (8*k*T./(pi.*mv)).^(1/2);
cp = (8*k*T/(pi*mp))^(1/2);

%mean free bath
lambda = 3.0.*(Diffp+D)./(cp.^2+cv.^2).^(1/2);  

%Diameter of vapor molecules
dv = (6.*mv./(pi.*rhoi)).^(1/3);

%Knudsen number
Kn = 2.*lambda./(dp+dv); 

%Transition regime correction factor for mass flux
beta = (1.0 + Kn)./(1.0 + (4./(3.0.*alpham) + 0.377).*Kn + 4.0./(3.0.*alpham).*Kn.^2.);

%List of all product compounds that are formed in this simulation
product_number=unique(product_numbers);

%--- Oligomerization and decomposition ---------------------------
%Molecular concentration in the particle (molecules/m^3)
CNpi = mpi./M*NA/Vp;
%Loss rates due to oligomerization (molecules/s)
L_olig = zeros(1,length(whats));
for ii=1:length(olig_ind)
    re_index=find(RxnP(:,1)==ii | RxnP(:,2)==ii);
    if RxnP(:,1)==ii & RxnP(:,2)==ii
    L_olig(olig_ind(ii)) = 2*nansum(RxnP(re_index,3).*(1-D_frac(olig_ind(RxnP(re_index,1)))).*...
       CNpi(olig_ind(RxnP(re_index,1))).*(1-D_frac(olig_ind(RxnP(re_index,2)))).*CNpi(olig_ind(RxnP(re_index,2))).*Vp);
    else
        L_olig(olig_ind(ii)) = nansum(RxnP(re_index,3).*(1-D_frac(olig_ind(RxnP(re_index,1)))).*...
       CNpi(olig_ind(RxnP(re_index,1))).*(1-D_frac(olig_ind(RxnP(re_index,2)))).*CNpi(olig_ind(RxnP(re_index,2))).*Vp);
    end
end

%Loss rates due to decompsition (molecules/s)
L_deco = zeros(1,length(whats));
if reversible==0
    for ii=1:length(decom_ind)
        re_index=find(DP(:,1)==ii);
        L_deco(decom_ind(ii)) = nansum(DP(re_index,2).*(1-D_frac(decom_ind(DP(re_index,1)))).*...
            CNpi(decom_ind(DP(re_index,1))).*Vp);
    end
else
    for ii=1:length(pro_ind)
        re_index=find(DP(:,1)==product_number(ii));
        L_deco(simulation_compounds==pro_ind(ii))=nansum(DP(re_index,2).*(1-D_frac(simulation_compounds==pro_ind(ii))).*...
            CNpi(simulation_compounds==pro_ind(ii)).*Vp);
    end
end
    
%Production rates due to oligomerization and decomposition (molecules/s)
P_olidec = zeros(1,length(whats));
reo_index=find(RxnP(:,4)==watermark | RxnP(:,5)==watermark);%find oligomerization reactions, that form water
    red_index=find(DP(:,3)==watermark | DP(:,4)==watermark);%find Decomposition reactions, that form water
   
       P_olidec(1) = nansum(RxnP(reo_index,3).*(1-D_frac(olig_ind(RxnP(reo_index,1)))).*...
       CNpi(olig_ind(RxnP(reo_index,1))).*(1-D_frac(olig_ind(RxnP(reo_index,2)))).*CNpi(olig_ind(RxnP(reo_index,2))).*Vp)...
       + nansum(DP(red_index,2).*(1-D_frac(decom_ind(DP(red_index,1)))).*CNpi(decom_ind(DP(red_index,1))).*Vp);
   
if reversible==0
    for ii=1:length(pro_ind)

        reo_index=find(RxnP(:,4)==product_number(ii) | RxnP(:,5)==product_number(ii));%find oligomerization reactions, that form this compound
        red_index=find(DP(:,3)==product_number(ii) | DP(:,4)==product_number(ii));%find Decomposition reactions, that form this compound
        if DP(:,3)==product_number(ii) & DP(:,4)==product_number(ii)
           P_olidec(length(whats)-length(pro_ind)+ii) = nansum(RxnP(reo_index,3).*(1-D_frac(olig_ind(RxnP(reo_index,1)))).*...
           CNpi(olig_ind(RxnP(reo_index,1))).*(1-D_frac(olig_ind(RxnP(reo_index,2)))).*CNpi(olig_ind(RxnP(reo_index,2))).*Vp)...
           + 2*nansum(DP(red_index,2).*(1-D_frac(decom_ind(DP(red_index,1)))).*CNpi(decom_ind(DP(red_index,1))).*Vp);
        else
            P_olidec(length(whats)-length(pro_ind)+ii) = nansum(RxnP(reo_index,3).*(1-D_frac(olig_ind(RxnP(reo_index,1)))).*...
           CNpi(olig_ind(RxnP(reo_index,1))).*(1-D_frac(olig_ind(RxnP(reo_index,2)))).*CNpi(olig_ind(RxnP(reo_index,2))).*Vp)...
           + nansum(DP(red_index,2).*(1-D_frac(decom_ind(DP(red_index,1)))).*CNpi(decom_ind(DP(red_index,1))).*Vp);
        end

    end
else
    for ii=1:length(pro_ind)

        reo_index=find(RxnP(:,4)==product_number(ii) | RxnP(:,5)==product_number(ii));%find oligomerization reactions, that form this compound
        red_index=find(DP(:,1)==product_number(ii));%find Decomposition reactions, that form this compound
        if DP(:,3)== DP(:,4)
           P_olidec(length(whats)-length(pro_ind)+ii) = nansum(RxnP(reo_index,3).*(1-D_frac(olig_ind(RxnP(reo_index,1)))).*...
           CNpi(olig_ind(RxnP(reo_index,1))).*(1-D_frac(olig_ind(RxnP(reo_index,2)))).*CNpi(olig_ind(RxnP(reo_index,2))).*Vp);
       P_olidec(olig_ind(DP(red_index,3)))=2*nansum(DP(red_index,2).*(1-D_frac(simulation_compounds==pro_ind(ii))).*CNpi(simulation_compounds==pro_ind(ii)).*Vp);
        else
            P_olidec(length(whats)-length(pro_ind)+ii) = nansum(RxnP(reo_index,3).*(1-D_frac(olig_ind(RxnP(reo_index,1)))).*...
           CNpi(olig_ind(RxnP(reo_index,1))).*(1-D_frac(olig_ind(RxnP(reo_index,2)))).*CNpi(olig_ind(RxnP(reo_index,2))).*Vp);
           P_olidec(olig_ind(DP(red_index,3:4)))= nansum(DP(red_index,2).*(1-D_frac(simulation_compounds==pro_ind(ii))).*CNpi(simulation_compounds==pro_ind(ii)).*Vp);
        end

    end
end

p_eq2=activ.*Xmole;
%-----------------------------------------------------------------

%--- Mass fluxes ----------------------------------------------------------
p_eq2(p_eq2<1e-15)=0;

flx_mass = M .* ( 2*pi*(dv+dp).*(D+Diffp).*beta.*(pvi-p_eq2)./(R*T) + (P_olidec'-L_olig'-L_deco')/NA );

flx_mass(mpi<=0 & (pvi-p_eq2)<0)=0;
% % Mass fluxes are calculated only for acids and neutrals (water is assumed to
% % be in equilibrium and moles of NH3 to equal moles of H2SO4)
flx = flx_mass(whats<0 | whats>=10);

%--------------------------------------------------------------------------

end
