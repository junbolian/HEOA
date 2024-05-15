%% Logistic chaotic mapping initialization
function Positions=initializationLogistic(pop,dim,ub,lb)

    Boundary_no= length(ub); % number of boundaries

    for i = 1:pop
        for j = 1:dim
            x0 = rand;
            a = 4;
            x = a*x0*(1-x0);
            if Boundary_no == 1
                Positions(i,j) = (ub-lb)*x + lb;
                if Positions(i,j)>ub
                    Positions(i,j) = ub;
                end
                if Positions(i,j)<lb
                    Positions(i,j) = lb;
                end
            else
                Positions(i,j) = (ub(j)-lb(j))*x + lb(j);
                if Positions(i,j)>ub(j)
                    Positions(i,j) = ub(j);
                end
                if Positions(i,j)<lb(j)
                    Positions(i,j) = lb(j);
                end
            end
            x0 = x;
        end
    end
end





