function[erpeeu, erperu, reponseInd, psi, tau] = performance(fct_transfert)
    %syms t
    % fonction de transfert comme argument
    % erpeeu : erreur en regime permanent entree echellon unite
    % erperu : erreur en regime permanent entree rampe pente unite
    % deno : polynome denominateur
    % numo : polynome numerateur
    % G(s) = [Nn Nn-1 Nn-2 ...N2 N1 N0 ; Dn Dn-1 Dn-2 ...D2 D1 D0]
    % reponseInd : fonction de t

    nbr_coeff = length(fct_transfert);
    deno = fct_transfert(2,1:nbr_coeff);
    numo = fct_transfert(1,1:nbr_coeff);
    D0= deno(end);
    D1 = deno(end-1);
    N0 = numo(end);
    type = 0;
    psi = Inf;

    if (D0 ~= 0)

        kp = N0/D0;
        erpeeu = 1 /(1+kp);
        erperu = Inf;
    elseif ((D0 ==0) && (D1 ~= 0))
        kv = N0/D1;
        erpeeu = 0;
        erperu = 1/kv;
        type = 1;
    elseif ((D0==0) && (D1==0))
        erpeeu = 0;
        erperu = 0;
        type = 2;
    end

    [reponseInd] = reponse(fct_transfert);
    tau = const(reponseInd)
end


%Fonction pour la reponse indicielle

function[fctRes] = reponse(fct_transfert)
    syms t
    %fct_transfert = fct_transfert./fct_transfert(2,1);
    N = cat(2, 0, fct_transfert(1,1:length(fct_transfert)))
    D = cat(2, fct_transfert(2,1:length(fct_transfert)), 0)
    [r, p, k] = residue(N, D)
    res = length(r);
    fctRes = 0;
    
    for i = 1:res
        fctRes = fctRes + r(i) * exp(p(i)*t) 
    end
end

% Lorsque nous avons des formes simples on utilise la fonction simple
function[reponseInd, psi] = simple(fct_transfert)
    syms t
    psi = Inf;
    if length(fct_transfert) == 2
        y = fct_transfert(1,2)/fct_transfert(2, 2) * (1 - exp(-fct_transfert(2, 2)*t));
        reponseInd = y;
    elseif length(fct_transfert) == 3
        gain = fct_transfert(1,3)/fct_transfert(2,end);
        fct_transfert(1, end) = fct_transfert(2,end);
        omega = sqrt(fct_transfert(1, end));
        psi = fct_transfert(2,end-1)/(omega * 2);
        theta = atan(psi/sqrt(1 - psi^2));
        y = gain*(1 - (exp(-psi*omega*t)/sqrt(1 - psi^2)) * sin(sqrt(1 - psi^2)*omega*t + theta));
        reponseInd = y;
    end
end

% Fonction pour caclculer les constantes de temps pour la question c
function[tmonte, tetabliss, dmax] = constantes(fct_transfert)
	taille = length(fct_transfert);
	tau = 1/fct_transfert(2,end);
	gain = fct_transfert(1,end);
	if (taille==2)
		tmonte = 2.2*tau;
		tetabliss = 4*tau;
		dmax = gain/fct_transfert(2,end);
	elseif (taille ==3)
		wncarre = fct_transfert(1,end);
		a = fct_transfert(2,2);
		amorti = a/sqrt(wncarre);
		
		tmonte = 1.8 / sqrt(wncarre);
		tetabliss = 8/a ;
		dmax = exp(-amorti*pi/sqrt(1-amorti^2))*100;
		
	end
end

function [tau] = const(fct)
    syms t;
    VF = limit(fct, t, Inf);
    VCT = 63/100*VF;
    type root2d.m;
    
    tau = fsolve(fct-VCT, [0 0]);
end



    

