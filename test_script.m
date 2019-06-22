%----------------------------------------------------------------------%
%---------------------------- Skript fuer -----------------------------%
%-------------- das modifizierte Polak-Ribiere-Verfahren --------------%
%----------- (am Beispiel der Rosenbrock- und Box-Funktion) -----------%
%------------------- im Vergleich zu PR, FR und HS --------------------%
%----------------------------------------------------------------------%

% Erstellen der Funktion und des Gradienten

fun = "box";                % Auswahl: "ros" fuer Rosenbrock, "box" fuer Box
switch fun
    case "ros"
        f = @ros;
        g = @gradros;
        x0 = [-1.2; 1];     % Wahl des Startpunktes
    case "box"
        f = @box;
        g = @gradbox;
        x0 = [0; 10; 20];   % Wahl des Startpunktes
end
    
% Bestimmung der Parameter
% Benutzer kann hier beliebige Werte im passenden Bereich waehlen

eps = 10^(-5);              % eps>=0
beta = 0.5;                 % 0<beta<1
sigma = 10^(-4);            % 0<sigma<1
delta1 = 0.1;               % 0<delta1<1
delta2 = 10;                % 1<delta2
maxsteps = 2000;            % Maximale Anzahl an Schritten

gamma = 0.1;                %0<gamma<1/2
eta = 0.7;                  %1/2<eta<1

% Modifiziertes Polak-Ribiere-Verfahren
[x_min, steps, counter] = modPRalg(maxsteps, x0, eps, beta, sigma, ...
                                                    delta1, delta2, f, g);
% Polak-Ribiere-Verfahren (mit PW-Regel)
[x_minPR, stepsPR, counterPR] = algFRHS(maxsteps, x0, eps, ...
                                                    eta, f, g, "PR");
% Fletcher-Reeves-Verfahren (mit PW-Regel)
[x_minFR, stepsFR, counterFR] = algFRHS(maxsteps, x0, eps, ...
                                                    eta, f, g, "FR");
% Hestenes-Stiefel-Verfahren (mit PW-Regel)
[x_minHS, stepsHS, counterHS] = algFRHS(maxsteps, x0, eps, ...
                                                    eta, f, g, "HS");
% % Gradientenverfahren (mit PW-Regel)
% [x_minGrad, stepsGrad, counterGrad] = algFRHS(maxsteps, x0, eps, ...
%                                                     eta, f, g, "Grad");


%----------------------------------------------------------------------%
%------------------------- Plot fuer Rosenbrock -----------------------%
%----------------------------------------------------------------------%

if (fun == "ros")
    
    % Gitterbildung
    x1 = linspace(-1.5, 1.5, 50);
    x2 = linspace(-0.5, 2, 50);
    [X1, X2]= meshgrid(x1, x2);
    
    % "Umdefinition" der Funktion, damit auf meshgrid anwendbar
    rostemp= 100*(X2 - X1.^2).^2 + (1 - X1).^2;
    
    figure
    
    % FR
    subplot(2,2,1)
    contour(X1, X2, rostemp, 75)                                           % plotten der Hoehenlinien
    hold on
    plot(stepsFR(:, 1), stepsFR(:, 2), 'o-g', 'LineWidth', 1.5)            % plotten der Iterierten (FR)
    plot(stepsFR(1, 1), stepsFR(1, 2), 'ok')                               % hervorheben: Startpunkt
    plot(stepsFR(end, 1), stepsFR(end, 2), 'xk')                           % hervorheben: Endpunkt

    % Beschriftung des Plots
    xlabel('x-Achse')
    ylabel('y-Achse')
    ylim([-0.5 2])
    title('Das Fletcher-Reeves-Verfahren')
    legend('Höhenlinien', 'x^k', 'x_0', '\approx (1,1)', 'Location', 'north')
    
    % PR
    subplot(2,2,2)
    contour(X1, X2, rostemp, 75)                                           % plotten der Hoehenlinie
    hold on
    plot(stepsPR(:, 1), stepsPR(:, 2), 'o-r', 'LineWidth', 1.5)            % plotten der Iterierten (PR)
    plot(stepsPR(1, 1), stepsPR(1, 2), 'ok')                               % hervorheben: Startpunkt
    plot(stepsPR(end, 1), stepsPR(end, 2), 'xk')                           % hervorheben: Endpunkt

    % Beschriftung des Plots
    xlabel('x-Achse')
    ylabel('y-Achse')
    title('Das Polak-Ribiere-Verfahren')
    legend('Höhenlinien', 'x^k', 'x_0', '\approx (1,1)', 'Location', 'north')
    
    
    % PR mod
    subplot(2,2,3)
    contour(X1, X2, rostemp, 75)                                           % plotten der Hoehenlinien
    hold on
    plot(steps(:, 1), steps(:, 2), 'o-b', 'LineWidth', 1.5)                % plotten der Iterierten (modPR)
    plot(steps(1, 1), steps(1, 2), 'ok')                                   % hervorheben: Startpunkt
    plot(steps(end, 1), steps(end, 2), 'xk')                               % hervorheben: Endpunkt

    % Beschriftung des Plots
    xlabel('x-Achse')
    ylabel('y-Achse')
    title('Das modifizierte Polak-Ribiere-Verfahren')
    legend('Höhenlinien', 'x^k', 'x_0', '\approx (1,1)', 'Location', 'north')
    
    
    % HS
    subplot(2,2,4)
    contour(X1, X2, rostemp, 75)                                           % plotten der Hoehenlinien
    hold on
    plot(stepsHS(:, 1), stepsHS(:, 2), 'o-m', 'LineWidth', 1.5)            % plotten der Iterierten (HS)
    plot(stepsHS(1, 1), stepsHS(1, 2), 'ok')                               % hervorheben: Startpunkt
    plot(stepsHS(end, 1), stepsHS(end, 2), 'xk')                           % hervorheben: Endpunkt

    % Beschriftung des Plots
    xlabel('x-Achse')
    ylabel('y-Achse')
    title('Das Hestenes-Stiefel-Verfahren')
    legend('Höhenlinien', 'x^k', 'x_0', '\approx (1,1)', 'Location', 'north')
    
    %disp(["FR", "PR", "modPR", "HS", "grad"; counterFR, counterPR, counter, counterHS])
    
    
    fprintf("#Schritte:\nFR\t PR\t modPR\t HS\t Grad (mit Armijo)\n%d\t %d\t %d\t %d\t %s\n", ...
                      counterFR, counterPR, counter, counterHS, ">10.000");
    
end