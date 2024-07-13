function runperf(fct_transfert, pas, limite)
    syms t
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
        [erpeeu, erperu, reponseInd, psi, tau, tr, te, depMax, p, vf, non_amortie] = performance(fct_transfert, pas, limite);
        if tau == Inf
            disp("Essayer de changer les valeurs du pas ou du temps d'études du système")
        elseif non_amortie == 0
            a = vf;
            reponseInd = simplify(reponseInd);
            fprintf("\nL'erreur en regime permanent pour une entrée échelon est : %e", erpeeu);
            fprintf("\nL'erreur en regime permanent pour une entrée rampe est : %e \n", erperu);
            disp('La reponse indicielle est :');
            disp(['y(t) = ', char(reponseInd)]);  % Convertir l'expression en chaîne de caractères

            
            fprintf("\nLa constant de temps est de : %e", tau);
            fprintf("\nLe temps de montée est de : %e sécondes", tr);
            fprintf("\nLe temps d'établissement est de : %e sécondes", te);
            
            if length(p) >= 3
                if round(depMax, 3) ~= 0 && length(p) >= 3
                    formatted_depMax = sprintf('%.1f', depMax);
                    fprintf("\nLe dépassement maximal en pourcentage est de  : ");
                    disp(formatted_depMax);
                    psi = -log(depMax/100)/(sqrt(pi^2 + (log(depMax/100))));
                elseif (depMax == 0 && length(p) == 3)
                    D = fct_transfert(2,1:length(fct_transfert));
                    D = D ./ D(1);
                    psi = D(2)/(2*sqrt(D(3)));
                end
            elseif length(p) <= 2
                disp('Le système est du premier ordre')
                psi = Inf;
            end
            courbe = ' SOUS-AMORTIE';
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
                
                fplot(reponseInd)
                %fplot(63/100)
                title(titre)
                grid on
                xlabel('temps t')
                ylabel('y(t)')
                axis([0 double(1.5*te) 0 double(a*2)])
                
            else
                courbe = 'PREMIER ORDRE';
                fplot(reponseInd)
                title 'Allure de la reponse indicielle du premier ordre'
                grid on
                xlabel('temps t')
                ylabel('y(t)')
                axis([0 double(1.5*te) 0 double(a*1.5)])
                
                
            end

            %% Génération du fichier Excel
            matrice= {'Title : ','Etudes des performances des systèmes automatiques';
                'Numerateur : ', num2str(fct_transfert(1,:));
                'Denominateur :', num2str(fct_transfert(2,:));
                'Le système est stable ? (1/0)', stable;
                'Nombre de pôles à gauche : ',gauche;
                'Nombre de pôles à gauche : ',droite;
                'Erreur en régime permanent pour une entrée échelon unitaire',erpeeu;
                'Erreur en régime permanent pour une entrée rampe de pente 1',erperu;
                'Constante de Temps',tau;
                'Nature du système : ',courbe;
                'Temps de montée : ',tr;
                "Temps d'établissement : ",te; 
                "Dépassement maximal en %: ", double(depMax);
                };
            writecell(matrice,'Rapport_Analyse.xls');
        
        else
                courbe = 'NON-AMORTIE';
                x_values = linspace(1, 5, 100);
                f(t) = reponseInd;
                y_values = double(f(x_values));
                fplot(reponseInd)
                titre = strcat('Réponse indicielle de nature : ', courbe);
                grid on
                title(titre)
                xlabel('temps t')
                ylabel('y(t)')
                ylim([min(y_values)*1.5 max(y_values)*1.5])
                 %% Génération du fichier Excel
                matrice= {'Title : ','Etudes des performances des systèmes automatiques';
                'Numerateur : ', num2str(fct_transfert(1,:));
                'Denominateur :', num2str(fct_transfert(2,:));
                'Le système est stable ? (1/0)', stable;
                'Nombre de pôles à gauche : ',gauche;
                'Nombre de pôles à gauche : ',droite;
                'Erreur en régime permanent pour une entrée échelon unitaire',erpeeu;
                'Erreur en régime permanent pour une entrée rampe de pente 1',erperu;
                'Constante de Temps', 'Non valide';
                'Nature du système : ', 'Non valide';
                'Temps de montée : ', 'Non valide';
                "Temps d'établissement : ",'Non valide'; 
                "Dépassement maximal : ", 'Non valide';
                };
            writecell(matrice,'Rapport_Analyse.xls');
                
        end
    end

    
end
