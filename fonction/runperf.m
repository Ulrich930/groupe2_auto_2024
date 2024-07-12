function runperf(fct_transfert)
    % fct_transfert point de départ du projet
    %% Critère de Routh
    % Tester la stabilité d'un système via sa fonction de transfer
    [stable, droite, gauche] = routh(fct_transfert(2,:));
    
    if stable == 0
        disp('Le système est instable selon le critère de Routh');
        fprintf('\nIl y a %g poles dans le demi-plan de gauche et %d dans le demi-plan de de droite\n', droite, gauche);
    else
        disp('Le sytème est stable selon le critère de Routh')
        disp('')
        disp('Nous pouvons passer à la prochaine étape')
        [erpeeu, erperu, reponseInd, psi, tau, tr, te, depMax, p] = performance(fct_transfert);
        reponseInd = simplify(reponseInd);
        fprintf("\nL'erreur en regime permanent pour une entrée échelon est : %e", erpeeu);
        fprintf("\nL'erreur en regime permanent pour une entrée rampe est : %e", erperu);
        disp('La reponse indicielle est :');
        disp(['y(t) = ', char(reponseInd)]);  % Convertir l'expression en chaîne de caractères

        
        fprintf("\nLa constant de temps est de : %e", tau);
        fprintf("\nLe temps de montée est de : %e sécondes", tr);
        fprintf("\nLe temps d'établissement est de : %e sécondes", te);
        if length(p) >= 3
            if depMax ~= 0 && length(p) >= 3
                fprintf("\nLe dépassement maximal est de : %e %", depMax);
                psi = -log(depMax/100)/(sqrt(pi^2 + (log(depMax/100))));
            elseif depMax == 0 && length(p) == 3
                D = cat(2, fct_transfert(2,1:length(fct_transfert)), 0);
                D = D ./ D(1);
                psi = D(2)/(2*D(3));
            end
        elseif length(p) <= 2
            fprinf('Le système est du premier ordre')
            psi = Inf;
        end
        if psi ~= Inf
            if psi == 0
                courbe = ' NON AMORTIE';
            elseif psi > 0 && psi < 1
                courbe = ' SOUS-AMORTIE';
            elseif psi == 1 
                courbe = ' CRITIQUE AMORTIE';
            elseif psi > 1
                courbe = ' SUR-AMORTIE';
            end
            titre = strcat('Réponse indicielle de nature : ', courbe);
            hold on
            fplot(reponseInd)
            %fplot(63/100)
            title(titre)
            grid on
            xlabel('temps t')
            ylabel('y(t)')
            axis([0 10 0 2])
            hold off
        else
            hold on
            fplot(reponseInd)
            title 'Allure de la reponse indicielle du premier ordre'
            grid on
            xlabel('temps t')
            ylabel('y(t)')
            hold off
        end
    end

end
