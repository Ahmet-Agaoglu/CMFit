function e_1 = UnknownErrorEstimate(x,lb,ub)
    % UnknownErrorEstimate: Returns the centroid of the FPS as the best estimate e_1.    

    % INPUT:
    % x  : Independent variable values
    % lb : Lower bounds of the dependent variable values
    % ub : Upper bounds of the dependent variable values

    % OUTPUT:
    % e_1 : centroid of the FPS [intercept, slope]

    % Find vertices of the FPS
    V = CorridorMethod(x,lb,ub);
    % Define the polygon (FPS) in the parameter space 
    S = polyshape(V(:,1),V(:,2));
    % Find the centroid of the polygon (used when distributional information is unavailable)
    [e_1(1), e_1(2)] = centroid(S);
    
    % Plot the result
    figure; hold on; plot(S); 
    plot(e_1(1),e_1(2),'o','MarkerSize', 8, 'MarkerFaceColor', 'r');
    legend('FPS','e_1')
    xlabel('intercept'); ylabel('slope'); grid on;
    hold off;
end
