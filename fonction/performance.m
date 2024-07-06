function[erpeeu, erperu] = performance(fct_transfert)
    % fonction de transfert comme argument
    % erpeeu : erreur en regime permanent entree echellon unite
    % erperu : erreur en regime permanent entree rampe unite
    % deno : polynome denominateur
    % T(s) = [Nn Nn-1 Nn-2 ...N2 N1 N0 ; Dn Dn-1 Dn-2 ...D2 D1 D0]

    deno = fct_transfert(2,:);
    numo = fct_transfert(:;2);
    D0= deno(end);
    D1 = deno(end-1);
    N0 = numo(end);
    N1 = numo(end-1)

    if D0 ~= 0
        kp = N0/D0
        erpeeu = 1 /(1+kp)
        erperu = Inf

    elseif (D0 ==0) && (D1 ~= 0)
        kv = N0/D1
        erpeeu = 0
        erperu = 1/kv

    elseif (D0==0) && (D1==0)
        erpeeu = 0
        erperu = 0

    end
end



    

