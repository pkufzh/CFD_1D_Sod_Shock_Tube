
% Function Module: Plot the properties with axis coordinate

function Plot_Props(t, xp, rho, p, u, E)

    subplot(2,2,1);
    plot(xp, rho, '-b', 'LineWidth', 2);
    xlabel('x');
    ylabel('Density');
    title([' T = ', num2str(t, '%.3f'), ' s']);
    axis([-0.5, 0.5, 0.0, 1.0]);
    grid on;
    
    subplot(2,2,2);
    plot(xp, p, '-g', 'LineWidth', 2);
    xlabel('x');
    ylabel('Pressure');
    title([' T = ', num2str(t, '%.3f'), ' s']);
    axis([-0.5, 0.5, 0.0, 1.0]);
    grid on;
    
    subplot(2,2,3);
    plot(xp, u, '-r', 'LineWidth', 2);
    xlabel('x');
    ylabel('Velocity');
    title([' T = ', num2str(t, '%.3f'), ' s']);
    axis([-0.5, 0.5, 0.0, 1.0]);
    grid on;
    
    subplot(2,2,4);
    plot(xp, E, '-m', 'LineWidth', 2);
    xlabel('x');
    ylabel('Specific Internal Energy');
    title([' T = ', num2str(t, '%.3f'), ' s']);
    axis([-0.5, 0.5, 1.5, 3.0]);
    grid on;

end

