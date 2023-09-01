%Calc_equil_water_AH 
function [p_eq,Masses]=Calc_equil_water(Masses)

global T R sigma P_sats pvi whats M psat_w V_molec NA
organs=find(whats==-2 | whats==-1 | whats==1 | whats==2 | whats>=10);
ite_dif=0.02;
psat_H2SO4=0;
psat_NH3=3.32e-6;
psat=zeros(size(whats));
psat(whats==0)=psat_w;
psat(whats==-3)=psat_H2SO4;
psat(whats==3)=psat_NH3;
psat(organs)=P_sats(organs);
Xmass=Masses/sum(Masses);
moles=Masses./M;
if any(whats==3)
    moles(whats==3)=moles(whats==-3);
end

activ=1;
  
 %Particle radius
 Vp=sum((moles)*NA.*V_molec);
 rp = (3*Vp/(4*pi)).^(1.0/3.0);

% Kelvin's effect
Ke = exp(2.*sigma.*V_molec.*NA/R./T./rp);  

still_not_close_enough=1;
while still_not_close_enough==1
    water_moles=(pvi(whats==0)*sum(moles(whats~=0)))/(activ*psat(whats==0)*Ke(whats==0)-pvi(whats==0));
    
    old_water_moles=moles(whats==0);
    moles(whats==0)=water_moles;
    
    Masses=moles.*M;
    Xmass=Masses/sum(Masses);
    
    %Particle radius
        Vp=sum((moles)*NA.*V_molec);
        rp = (3*Vp/(4*pi)).^(1.0/3.0);

    % Kelvin's effect
    Ke = exp(2.*sigma.*V_molec.*NA/R./T./rp);  
    
    water_dif=abs((water_moles-old_water_moles)./water_moles);
    
    if water_dif<=ite_dif
        still_not_close_enough=0;
    end
end

Xmol=moles/sum(moles([find(whats==0); organs]));
p_eq=(activ.*Xmol.*psat.*Ke);
 






