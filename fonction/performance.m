function[erpeeu, erperu, reponseInd, psi, tau, tr, te, depMax, p, vf, non_amortie] = performance(fct_transfert, pas, limite)
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
    [tau, tr, te, depMax, vf, non_amortie] = const(reponseInd, pas, limite);
end


%Fonction pour la reponse indicielle

function[fctRes, p] = reponse(fct_transfert)
    syms t;
    %fct_transfert = fct_transfert./fct_transfert(2,1);
    N = cat(2, 0, fct_transfert(1,1:length(fct_transfert)));
    D = cat(2, fct_transfert(2,1:length(fct_transfert)), 0);
    [r, p] = residue(N, D);

    polUt = [];
    multi = 1;
    res = length(r);
    fctRes = 0;
    for i = 1:res
        if ismember(p(i), polUt) 
            fctRes = fctRes + r(i)/(factorial(multi))*t^(multi)*exp(p(i)*t);
            multi = multi + 1;
        else
            multi = 1;
            fctRes = fctRes + r(i) * exp(p(i)*t);
            if imag(p(i)) == 0
                polUt = [polUt, [p(i)]];
            end
        end
        
    end
end

% Lorsque nous avons des formes simples on utilise la fonction simple
function[reponseInd, psi, omega] = simple(fct_transfert)
    omega = 0;
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

function [tau, tr, te, depMax, vf, non_amortie] = const(fct, pas, limite)
    try
        syms t;
            non_amortie = 0;
            % Tentative d'exécution d'une opération risquant de générer une erreur
            vf = limit(fct, t, Inf);
            f(t) = fct;
            fplot(fct, [0 10^-2]);
        if ~isnan(vf)
            % Calcule de la limite de fct lorsque t tend vers l’infini, stockée dans vf
            if vf >= 1
                VCT = round(63/100*vf, 2);
                VTM1 = round(0.1*vf, 2);
                VTM2 = round(0.9*vf, 2);
                x_values = linspace(0, pas, limite);
                y_values = round(double(f(x_values)), 2);
                a = find(y_values == VCT);
                b1 = find(y_values == VTM1);
                b2 = find(y_values == VTM2);
                tr1 = x_values(b1(1));
                tr2 = x_values(b2(1));
                tr = tr2 - tr1;
                depMax = simplifyFraction(((max(y_values) - vf)/vf)*100);
                VTE1 = round(vf*98/100, 2);
                VTE2 = round(vf+ vf*0.02, 2);
                c1 = find(round(y_values, 2) == VTE1);
                c2 = find(round(y_values, 2) == VTE2);
                if isempty(c2)
                    te = x_values(c1(end));
                else
                    te = max(x_values(c1(end)), x_values(c2(end)));
                end
            else
                i = 0;
                while round(vf, i) == 0
                    i = i + 1;
                end
                VCT = round(63/100*vf, i+1);
                VTM1 = round(0.1*vf, i+1);
                VTM2 = round(0.9*vf, i+1);
                x_values = linspace(0, pas, limite);
                y_values = double(f(x_values));
                round(y_values(1:100), i+1)
                a = find(round(y_values, i+1) == VCT);
                searchIter = 1;
                while isempty(a) 
                    VCT = VCT-0.1^(i+1);
                    a = find(round(y_values, i+1) == VCT); 
                    searchIter  = searchIter + 1;
                    if searchIter == 5
                        disp("Modifier la valeur du pas ou de la limite du temps d'étude")
                        break
                    end
    
                end
                b1 = find(round(y_values, i+1) == VTM1);
                b2 = find(round(y_values, i+1) == VTM2);
                while isempty(b2)
                    VTM2 = VTM2-0.1^(i+1);
                    b2 = find(round(y_values, i+1) == VTM2);
                end
                if isempty(b1)
                    tr1 = 0;
                    tr2 = x_values(b2(1));
                else
                    tr1 = x_values(b1(1));
                    tr2 = x_values(b2(1));
                end
                tr = tr2 - tr1;
                depMax = simplifyFraction(((max(y_values) - vf)/vf)*100);
                if depMax < vf
                    depMax =0;
                end
                VTE1 = round(vf*98/100, i+1);
                VTE2 = round(vf+ vf*0.02, i+1);
                c1 = find(round(y_values, i+1) == VTE1);
                while isempty(c1)
                    VTE1 = VTE1-0.1^(i+1);
                    c1 = find(round(y_values, i+1) == VTE1);  
                end
                c2 = find(round(y_values, i+1) == VTE2);
                if isempty(c2)
                    te = x_values(c1(end));
                else
                    te = max(x_values(c1(end)), x_values(c2(end)));
                end
            end
            tau = x_values(a(1));
        else
            disp('Ici')
            x_values = linspace(0, pas, limite);
            y_values = round(double(f(x_values)), 2);
            if length(find(y_values == max(y_values))) > 1
                non_amortie = 1;
                tau = 0;
                tr = 0;
                te = 0;
                depMax = 0;
                vf = 0;
            end
                
        end
   catch
            tau = Inf;
            tr = 0;
            te = 0;
            depMax = 0;
            vf = 0;  
    end
end



    

