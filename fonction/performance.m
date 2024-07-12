function[erpeeu, erperu, reponseInd, psi, tau, tr, te, depMax, p] = performance(fct_transfert)
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
    psi = Inf;
    if (D0 ~= 0)
        kp = N0/D0;
        erpeeu = 1 /(1+kp);
        erperu = Inf;
    elseif ((D0 ==0) && (D1 ~= 0))
        kv = N0/D1;
        erpeeu = 0;
        erperu = 1/kv;
    elseif ((D0==0) && (D1==0))
        erpeeu = 0;
        erperu = 0;
    end

    [reponseInd, p] = reponse(fct_transfert);
    [tau, tr, te, depMax] = const(reponseInd);
end


%Fonction pour la reponse indicielle

function[fctRes, p] = reponse(fct_transfert)
    syms t
    %fct_transfert = fct_transfert./fct_transfert(2,1);
    N = cat(2, 0, fct_transfert(1,1:length(fct_transfert)));
    D = cat(2, fct_transfert(2,1:length(fct_transfert)), 0);
    [r, p, k] = residue(N, D);
    res = length(r);
    fctRes = 0;
    
    for i = 1:res
        fctRes = fctRes + r(i) * exp(p(i)*t);
    end
end

% Lorsque nous avons des formes simples on utilise la fonction simple
function[reponseInd, psi] = simple(fct_transfert)
    % fonction nommée 'simple' qui prend en entrée 'fct_transfert'
    % et retourne 'reponseInd' et 'psi'.
    syms t
    psi = Inf;
    % Initialisation de  psi à l’infini.
    if length(fct_transfert) == 2
        y = fct_transfert(1,2)/fct_transfert(2, 2) * (1 - exp(-fct_transfert(2, 2)*t));
        % Calcule de  y en utilisant les éléments de fct_transfert et une fonction exponentielle.
        reponseInd = y;
    elseif length(fct_transfert) == 3
        gain = fct_transfert(1,3)/fct_transfert(2,end);
        % Calcule le gain en utilisant les éléments de fct_transfert
        fct_transfert(1, end) = fct_transfert(2,end);
        omega = sqrt(fct_transfert(1, end));
        psi = fct_transfert(2,end-1)/(omega * 2);
        theta = atan(psi/sqrt(1 - psi^2));
        y = gain*(1 - (exp(-psi*omega*t)/sqrt(1 - psi^2)) * sin(sqrt(1 - psi^2)*omega*t + theta));
        reponseInd = y;
    end
end

function [tau, tr, te, depMax] = const(fct)
    syms t;
    vf = limit(fct, t, Inf);
    % Calcule de la limite de fct lorsque t tend vers l’infini, stockée dans vf
    VCT = 63/100*vf;
    VTM1 = 0.1*vf;
    VTM2 = 0.9*vf;
    VTE1 = vf*98/100;
    VTE2 = vf+ vf*0.02;
    f(t) = fct;
    x_values = linspace(0, 5, 5000);
    y_values = double(f(x_values));
    x = round(x_values, 2);
    y = round(y_values, 2);
    
    a = find(y == VCT);
    b1 = find(y == VTM1);
    b2 = find(y == VTM2);
    c1 = find(y == VTE1);
    c2 = find(y == VTE2);
    tr1 = x_values(b1(round(length(b1)/2)));
    tr2 = x_values(b2(round(length(b2)/2)));
    tr = tr2 - tr1;
    depMax = simplifyFraction(((max(y_values) - vf)/vf)*100);
    te = max(x(c1(end)), x(c2(end)));
    tau = x_values(a(round(length(a)/2)));
end



    

