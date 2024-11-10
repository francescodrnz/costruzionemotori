function [delta_p] = deltaPFun(D_A,D_fori,dk,l_k,E_p)
%deltaPFun calcolo cedevolezza vite
tan_phi=0.362+0.032*log(l_k/(2*dk))+0.153*log(D_A/dk);

D_Alim=dk+l_k*tan_phi; % w = 1 sistema vite+dado

% calcolo delta_p

if D_A>=D_Alim
    
    delta_p=2/(pi*E_p*tan_phi*D_fori)*log(((dk+D_fori)*(dk+l_k*tan_phi-D_fori))/((dk-D_fori)*(dk+l_k*tan_phi+D_fori)));
    
elseif D_A<D_Alim && D_A>=dk
    
    delta_p=1/(pi*E_p)*(2/(D_fori*tan_phi)*log(((dk+D_fori)*(D_A-D_fori))/((dk-D_fori)*(D_A+D_fori)))+4/(D_A^2-D_fori^2)*(l_k-(D_A-dk)/(tan_phi)));
    
elseif D_A<dk
    
    A_p=pi/4*(dk^2-D_A^2);
    delta_p=l_k/(E_p*A_p);
   
end

end