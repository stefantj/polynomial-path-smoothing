using PyPlot # For plotting

type PolyProblem
    # Constraints
    B_x::Vector{Float64}        # Vector of constraint values for x dimension
    B_y::Vector{Float64}        # Vector of constraint values for y dimension
    B_orders::Vector{Int64}     # Vector of orders of constraints
    B_time_inds::Vector{Int64}  # Vector of indices of time for constraints
    # Cost behavior
    q_coeff::Vector{Float64}    # Cost vector
    kT::Float64                 # Time cost
end

type Poly2Segs                  # For 1 segment per pair of points
    prob::PolyProblem           # Original poly problem
    num_segs::Int64             # Number of segments
    cont_order::Int64           # Continuity order between segments
    degrees::Vector{Int64}      # Degree for each segment
    M::Array{Float64}           # Mapping matrix
end


# Forms a vanilla problem
function random_poly_prob(num_points::Int64)
    # Random points
    xpts = num_points*rand(num_points);#(0:num_points-1);
    ypts = num_points*rand(num_points);
    # Initial derivative constraints: 
    #          vel acc jerk
    bx_init = [1.0; 0.0; 0.0];
    by_init = [0.0; 0.0; 0.0];
    # Final derivative constraints:
    bx_final = [1.0; 0.0; 0.0];
    by_final = [0.0; 0.0; 0.0];

    #### Assemble constraint variables: #### 
    # Total B vectors:
    B_x = [xpts; bx_init; bx_final];
    B_y = [ypts; by_init; by_final];
    # Order of constraints, in order of B_x
    B_orders = [zeros(num_points); 1;2;3;1;2;3];
    # Time indices of constraints:
    B_time_inds = [collect(1:num_points); ones(3); num_points*ones(3)];

    degree = num_points+3+3;
    #### Define our cost metrics ####
    # Q coeffs: weight which derivatives you care about:
    q_coeff = zeros(degree);
    q_coeff[4] = 1;
    # Total time penalty
    kT = 5000; 
    
    return PolyProblem(B_x,B_y,round(Int64,B_orders),round(Int64,B_time_inds),q_coeff,kT)
end

# Forms segments between each pair of points in PolyProblem
function form_2segs(prob::PolyProblem)
    # Number of segments is number of points - 1:
    num_pts = size(find(prob.B_orders .== 0),1);
    if(num_pts != maximum(prob.B_time_inds))
        println("Error in counting number of points!")
    end

    num_segs = num_pts-1;

    cont_order = 5

    # Compute number of independent coefficients for each segment:
    # Most have degree cont_order
    indep_degrees=cont_order*ones(num_segs);
    # First poly has order cont_order + init_order
    indep_degrees[1] =  cont_order + size(find(prob.B_time_inds.== 1),1); 
    # Last poly has order final_order:
    indep_degrees[end] = size(find(prob.B_time_inds.== num_pts),1)
    # All other degrees are cont_order
    

end



function print_poly_prob(p::PolyProblem)
    println("********************************");
    println("PolyProblem: ")
    println("Constraints:")
    println("Order\tTime\tX\tY");
    for k=1:size(p.B_x,1)
        println(p.B_orders[k],"\t",p.B_time_inds[k],"\t",round(p.B_x[k],2),"\t",round(p.B_y[k],2));
    end
    println("Time penalty: $kT")
    println("Q_coeffs: $q_coeff");
    println("********************************");
end


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
        n=round(Int64,n);
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
        A = [A; constr_order(orders[k], time_vec[round(Int64,time_inds[k])], degree)];
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
    return J[1],x_coeff,y_coeff;
end


function solve_poly_problem(prob::PolyProblem, times::Vector{Float64})


    # Form A matrix:
    A = zeros(0,prob.degree);
    for seg=1:prob.num_segs
        for k=1:size(prob.B_x,1)
            A = [A; constr_order(prob.B_orders[k],times[prob.B_time_inds[k]],prob.degrees[seg])];
        end
    end

    # Need to form M matrix:


    # Form Q matrix:
    Q_mat = form_Q(prob.B_q_coeff, time_vec[end])

    # Compute coefficients:
    x_coeff = A\prob.B_x;
    y_coeff = A\prob.B_y;

    # Calculate objective
    J = x_coeff'*Q_mat*x_coeff + y_coeff'*Q_mat*y_coeff + kT*time_vec[end];

    return J, x_coeff, y_coeff
end
