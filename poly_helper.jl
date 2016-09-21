using PyPlot # For plotting

# Function for a prettier plot of polynomials with first few derivatives
function plot_poly(cx,cy,pp, T0, TF, fn)
    figure(fn);
    radius = 50;
    degree = length(cx);
    times = linspace(T0,TF,100)
    x0 = zeros(length(times))
    y0 = zeros(length(times))
    x1 = zeros(length(times))
    y1 = zeros(length(times))
    x2 = zeros(length(times))
    y2 = zeros(length(times))
    x3 = zeros(length(times))
    y3 = zeros(length(times))
    k = 0
    for t in times
        k+=1;
        for d = 1:degree
            x0[k] += cx[d] * t^(d-1)
            y0[k] += cy[d] * t^(d-1)
            if(d > 1)
                x1[k] += (d-1)*cx[d] * t^(d-2)
                y1[k] += (d-1)*cy[d] * t^(d-2)
            end
            if(d > 2)
                x2[k] += (d-1)*(d-2)*cx[d] * t^(d-3)
                y2[k] += (d-1)*(d-2)*cy[d] * t^(d-3)
            end
            if(d > 3)
                x3[k] += (d-2)*(d-1)*(d-3)*cx[d] * t^(d-4)
                y3[k] += (d-2)*(d-1)*(d-3)*cy[d] * t^(d-4)
            end
        end
    end

    subplot(2,2,1); title("Position");
    if(T0 == 0)
        plot(x0,y0, color=:green);
    else
        plot(x0,y0,color=:red);
    end
#    axis([-5,2*radius+5,-5,2*radius+5])
#    axis([-5,65,-5,35])
    scatter(pp[:,1], pp[:,2]);

    subplot(2,2,2); title("Velocity")
    plot(times, x1, color=:red)
    plot(times, y1, color=:blue)
    legend(["X", "Y"]);

    subplot(2,2,3); title("Acceleration")
    plot(times, x2, color=:red)
    plot(times, y2, color=:blue)
    legend(["X", "Y"]);

    subplot(2,2,4); title("Jerk");
    plot(times, x3, color=:red)
    plot(times, y3, color=:blue)
    legend(["X", "Y"]);
end



# Function for making Jacobian for free variable optimization
function form_Q(Q_coeffs, t)
    degree = size(Q_coeffs,1)
    Q_mat = zeros(degree,degree)
    for k = 0:degree-1
        if(Q_coeffs[k+1] == 0)
            continue;
        end
        c_k = Q_coeffs[k+1];
        for i=1:degree
        # Form Q
            for l = 1:degree
                if(i >= k && l >= k)
                    c_tmp = 2*c_k
                    for m = 0:k-1
                        c_tmp = c_tmp * (i-m)*(l-m);
                    end
                    Q_mat[i,l] += c_tmp*(t^(i+l-2*k+1))/(i+l-2*k+1)
                end
            # Else 0
            end
        end
    end
    return Q_mat
end

# Function for deriving the row of A corresponding to constraint of order at time for poly of degree.
function constr_order( order, time, degree)
    if(order==0) 
        row_vec = ones(1,degree); 
        for n=2:degree
            row_vec[1,n] = time^(n-1);
        end
        return row_vec;
    end
    row_vec = zeros(1,degree);
    for n = order:degree-1
        coeff = 1; 
        for k=1:order
            coeff = coeff*(n+1-k)
        end

        row_vec[1,n+1] = coeff * (time^(n-order))
    end
    return row_vec;
end


function solve_poly_problem( orders, time_inds, B_x, B_y, time_vec, q_coeffs,k_T)
## Form A matrix:
    A = zeros(0,degree)
    # Point constraints:    
    for k=1:size(B_x,1)
        A = [A; constr_order(orders[k], time_vec[time_inds[k]], degree)];
    end
    # Form Q matrix:
    Q_mat = form_Q(q_coeff, time_vec[end]);
    
    # for now, ignore explicit snap minimization
    # Otherwise should compute -Rpp\Rfp and find the B_free vector

    
    # Solve for coefficients:
    x_coeff = A\B_x;
    y_coeff = A\B_y;
    
    # Evaluate cost:
    J = x_coeff'*Q_mat*x_coeff + y_coeff'*Q_mat*y_coeff + kT*time_vec[end];
    return J[1];
end



