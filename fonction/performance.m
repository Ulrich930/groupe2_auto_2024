function[erpeeu, erperu, reponseInd, psi] = performance(fct_transfert)
    %syms t
    % fonction de transfert comme argument
    % erpeeu : erreur en regime permanent entree echellon unite
    % erperu : erreur en regime permanent entree rampe pente unite
    % deno : polynome denominateur
    % numo : polynome numerateur
    % T(s) = [Nn Nn-1 Nn-2 ...N2 N1 N0 ; Dn Dn-1 Dn-2 ...D2 D1 D0]
    % reponseInd : fonction de t

    nbr_coeff = length(fct_transfert);
    deno = fct_transfert(2,1:nbr_coeff);
    numo = fct_transfert(1,1:nbr_coeff);
    D0= deno(end);
    D1 = deno(end-1);
    N0 = numo(end);
    type = 0;

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

    [reponseInd, psi] = reponse(fct_transfert, type);
end

function[reponseInd, psi] = reponse(fct_transfert, type)
    psi = Inf;
    syms t
    fct_transfert = fct_transfert./fct_transfert(2,1);
    if type == 0
        if length(fct_transfert) == 2
            y = fct_transfert(1,2)/fct_transfert(2, 2) * (1 - exp(-fct_transfert(2, 2)*t));
            reponseInd = y;
        elseif length(fct_transfert) == 3
            gain = fct_transfert(1,3)/fct_transfert(2,end)
            fct_transfert(1, end) = fct_transfert(2,end);
            omega = sqrt(fct_transfert(1, end));
            psi = fct_transfert(2,end-1)/(omega * 2);
            theta = atan(psi/sqrt(1 - psi^2));
            y = gain*(1 - (exp(-psi*omega*t)/sqrt(1 - psi^2)) * sin(sqrt(1 - psi^2)*omega*t + theta));
            reponseInd = y;
        end
    end
end



    

