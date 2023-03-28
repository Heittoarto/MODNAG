%Calc_equil_water_AH 
function [p_eq,Masses,err,Xmol,Ke]=Calc_equil_water_AH_Bolivia(Masses)

global T R sigma P_sats pvi whats M psat_w rhoi
organs=find(whats==-2 | whats==-1 | whats==1 | whats==2 | whats>=10);
ite_dif=0.02;
psat_H2SO4=0;
psat_NH3=3.32e-6;
psat=zeros(size(whats));
psat(whats==0)=psat_w;
psat(whats==-3)=psat_H2SO4;
psat(whats==3)=psat_NH3;
psat(organs)=P_sats*101325;
Xmass=Masses/sum(Masses);
moles=Masses./M;
if any(whats==3)
    moles(whats==3)=moles(whats==-3);
end

activ=1;
  
 mp=sum(Masses);
% Particle density
rho = ( sum(Xmass./rhoi) ).^(-1);
       
% Particle radius        
rp = (3*(mp/rho)/(4*pi)).^(1.0/3.0);

% Kelvin's effect
Ke = exp(2.*M.*sigma./R./T./rp./rho);

still_not_close_enough=1;
loop_nr=0;
err=0;
while still_not_close_enough==1
    loop_nr=loop_nr+1;
    water_moles=(pvi(whats==0)*sum(moles(whats~=0)))/(activ*psat(whats==0)*Ke(whats==0)-pvi(whats==0));
    
    old_water_moles=moles(whats==0);
    moles(whats==0)=water_moles;
    
    Masses=moles.*M;
    Xmass=Masses/sum(Masses);
    
    % Particle density
    rho = ( sum(Xmass./rhoi) )^(-1);
       
    mp=sum(Masses);
    % Particle radius        
    rp = (3*(mp/rho)/(4*pi)).^(1.0/3.0);

    % Kelvin's effect
	Ke = exp(2.*M.*sigma./R./T./rp./rho);
    
    water_dif=abs((water_moles-old_water_moles)./water_moles);
    
    if water_dif<=ite_dif
        still_not_close_enough=0;
    end
    
    if loop_nr==15
        warning('myfun:warncode','Allowed iteration difference is too small for code to move forward') 
        ite_dif=ite_dif+0.01;
        err=1;
    elseif loop_nr==23
        ite_dif=ite_dif+0.01;
        err=2;
    elseif loop_nr==31
        ite_dif=ite_dif+0.01;
        err=3;
    elseif loop_nr==39
        ite_dif=ite_dif+0.01;
        err=4;
    elseif loop_nr==47
        ite_dif=ite_dif+0.01;
        err=5;
    elseif loop_nr==55
        ite_dif=ite_dif+0.01;
        err=6;
    elseif loop_nr==63               
        err=7;
    break
    end 
end


% Xmol=moles/sum(moles([find(whats==0); organs]));
Xmol=moles/sum(moles);
% Masses=Masses;
p_eq=(activ.*Xmol.*psat.*Ke);
end
 






