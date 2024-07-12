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
        [erpeeu, erperu, reponseInd, psi] = performance(fct_transfert);
        reponseInd = simplify(reponseInd);
        fprintf("\nL'erreur en regime permanent pour une entrée échelon est : %e", erpeeu);
        fprintf("\nL'erreur en regime permanent pour une entrée rampe est : %e", erperu);
        disp('La reponse indicielle est :');
        disp(['y(t) = ', char(reponseInd)]);  % Convertir l'expression en chaîne de caractères

        if psi == Inf
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
            title(titre)
            grid on
            xlabel('temps t')
            ylabel('y(t)')
            axis([0 10 0 2])
            hold off
        else
            hold on
            fplot(reponseInd)
            title 'Allure de la reponse indicielle'
            grid on
            xlabel('temps t')
            ylabel('y(t)')
            hold off
        end
        
         
    end


