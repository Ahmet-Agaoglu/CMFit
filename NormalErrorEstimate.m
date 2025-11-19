function e_2 = NormalErrorEstimate(x,y,lb,ub)
    % NormalErrorEstimate: Returns the optimal solution e_2 that maximizes likelihood function within the FPS.

    % INPUT:
    % x  : Independent variable values
    % y  : Response variable measurements
    % lb : Lower bounds of the dependent variable values
    % ub : Upper bounds of the dependent variable values

    % OUTPUT:
    % e_2 : optimal solution [intercept, slope]


    V = CorridorMethod(x,lb,ub); % Get the vertices of the FPS
    C = [];
    OLS = polyfit(x,y,1); % Calculate OLS solution
    OLS_Vals = OLS(1)*x+OLS(2); % Evaluate OLS solution
    if OLS_Vals >= lb & OLS_Vals <= ub
        % If OLS lies within FPS, it is the optimal solution
        e_2(1)=OLS(2); e_2(2)=OLS(1);
    else
        % Iterate over each edge of the polygon
        for i = 1:size(V,1)
            if i<size(V,1)
                V1 = V(i,:); V2 = V(i+1,:);
            else
                V1 = V(i,:); V2 = V(1,:);
            end
            % Minimize objective along the edge
            p = (V2(2)-V1(2))/(V2(1)-V1(1));
            q = V1(2)-p*V1(1);
            B0_bar = sum((y-q*x).*(1+p*x))/sum((1+p*x).^2);
            B1_bar = sum(p*B0_bar+q);
    
            % Check if the point lies within the edge
            if B0_bar<=max(V1(1),V2(1)) && B0_bar>=min(V1(1),V2(1)) && B1_bar<=max(V1(2),V2(2)) && B1_bar>=min(V1(2),V2(2))
                C = [C;B0_bar B1_bar sum((y-B0_bar-B1_bar*x).^2)];
            else
                % Evaluate objective at edge endpoints
                obj1 = sum((y-V1(1)-V1(2)*x).^2);
                obj2 = sum((y-V2(1)-V2(2)*x).^2);
                if obj1<obj2
                    C = [C; V1 obj1];
                else
                    C = [C; V2 obj2];
                end
            end
        end
        % Find the optimal point with the minimum objective value
        C = sortrows(C,3,"ascend");
        e_2(1)= C(1,1); e_2(2) = C(1,2);
    end
    % Plot the result
    % Define the polygon (FPS) in the parameter space 
    S = polyshape(V(:,1),V(:,2));
    figure; hold on; plot(S); 
    plot(e_2(1),e_2(2),'o','MarkerSize', 8, 'MarkerFaceColor', 'r');
    legend('FPS','e_2')
    xlabel('intercept'); ylabel('slope'); grid on;
    hold off;
end
