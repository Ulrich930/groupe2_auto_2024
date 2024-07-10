function runperf(fct_transfert)
    % fct_transfert point de départ du projet
    %% Critère de Routh
    % Tester la stabilité d'un système via sa fonction de transfer
    [stable, droite, gauche] = routh(fct_transfert(2,:));
    
    if stable == 0
        disp('\nLe système est instable selon le critère de Routh');
        fprintf('\nIl y a %g poles dans le demi-plan de gauche et %d dans le demi-plan de de droite\n', droite, gauche);
    else
        disp('Le sytème est stable selon le critère de Routh')
        disp('')
        disp('Nous pouvons passer à la prochaine étape')
    end


