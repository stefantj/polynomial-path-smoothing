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
    poly_prob::PolyProblem      # Original poly problem
    times::Vector{Float64}      # Times for each segment
    tot_degree::Int64           # Total reduced degree 
    num_segs::Int64             # Number of segments
    cont_order::Int64           # Continuity order between segments
    degrees::Vector{Int64}      # Degree for each segment
    C_mats::Array{Float64,3}    # Vector of matrices for converting segments polynomials into reduced polynomials
    x_coeffs::Array{Float64,2}  # X coefficients
    y_coeffs::Array{Float64,2}  # Y coefficients
end


# Forms a vanilla problem
function random_poly_prob(num_points::Int64)
    # Random points
    xpts = num_points*rand(num_points);#(0:num_points-1);
    ypts = num_points*rand(num_points);

#    xpts = collect(1:num_points);
#    ypts = num_points*sin(pi*(1:num_points)/5)

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
    q_coeff[1] = 1.0;
    q_coeff[2] = 0.5;
    q_coeff[3] = 0.1;
    q_coeff[4] = 1.0;
    # Total time penalty
    kT = 50000; 
    
    return PolyProblem(B_x,B_y,round(Int64,B_orders),round(Int64,B_time_inds),q_coeff,kT)
end

# Forms segments between each pair of points in PolyProblem
function form_2segs(prob::PolyProblem,times)
    # Number of segments is number of points - 1:
    num_pts = size(find(prob.B_orders .== 0),1);
    if(num_pts != maximum(prob.B_time_inds))
        println("Error in counting number of points!")
    end
    num_init_constr = size(find(prob.B_time_inds.==1),1);
    num_fin_constr  = size(find(prob.B_time_inds.==num_pts),1);
    num_segs = round(Int64,num_pts-1);

    cont_order = 5

    # Compute number of independent coefficients for each segment:
    # Most have degree 1 (for their 1 point)
    indep_degrees=round(Int64,ones(num_segs))
    # First poly has order 1 + init_order
    indep_degrees[1] =  round(Int64, num_init_constr + 1);
    # Last poly has order final_order:
    indep_degrees[end] = round(Int64, num_fin_constr);

    # Compute number of dependent coefficients for each segment:
    # Most have cont_order dependent coeffs
    dep_degrees = round(Int64,cont_order*ones(num_segs));
    # Initial is the one exception, and has no dependent coefficients.
    dep_degrees[1] = 0;

    degrees = indep_degrees+dep_degrees;

    # Size of the reduced polynomial coefficients
    tot_degree = round(Int64,sum(indep_degrees))
    # Here we compute the M recursively:
    C_mats = zeros(tot_degree,tot_degree,num_segs);
    C_mats[:,:,1] = [[eye(indep_degrees[1]) zeros(indep_degrees[1],tot_degree-indep_degrees[1])]; zeros(tot_degree-indep_degrees[1], tot_degree)];

    # Backsolve for C_mat[:,:,seg]
    for seg=2:num_segs
        indep_before = sum(indep_degrees[1:seg-1]);
        indep_after = sum(indep_degrees[seg+1:num_segs])
        Cprev = C_mats[:,:,seg-1];
        Cseg = [zeros(cont_order, tot_degree); [zeros(indep_degrees[seg],indep_before) eye(indep_degrees[seg]) zeros(indep_degrees[seg],indep_after)]];


        totdeg_before = indep_degrees[seg-1]+dep_degrees[seg-1];
        totdeg_here = indep_degrees[seg]+dep_degrees[seg];

        max_deg = max(totdeg_before-1,totdeg_here-1)
        for k=cont_order:-1:1
            Cseg[k,:] = Cprev[k,:];
            t_init = times[seg];
            for i=k:max_deg
                coeff = t_init^(i-k+1)
                for m=1:(k-1)
                    coeff*= (i-(m-1))/m;
                end
                if(i < totdeg_before)
                    Cseg[k,:] += coeff*Cprev[i+1,:];
                end
                if(i < totdeg_here)
                    Cseg[k,:] -= coeff*Cseg[i+1,:];
                end
            end
        end 
        println("Condition number for segment $seg: ", cond(Cseg));
        C_mats[1:totdeg_here,:,seg] = Cseg;
    end

    x_coeffs = zeros(tot_degree,num_segs)
    y_coeffs = zeros(tot_degree,num_segs)

   return Poly2Segs(prob, times, tot_degree, num_segs, cont_order, degrees, C_mats, x_coeffs, y_coeffs); 
end


function print_poly_prob(p::PolyProblem)
    println("********************************");
    println("PolyProblem: ")
    println("Constraints:")
    println("Order\tTime\tX\tY");
    for k=1:size(p.B_x,1)
        println(p.B_orders[k],"\t",p.B_time_inds[k],"\t",round(p.B_x[k],2),"\t",round(p.B_y[k],2));
    end
    println("Time penalty: ",p.kT)
    println("Q_coeffs: ",p.q_coeff);
    println("********************************");
end


function plot_poly(prob::Poly2Segs, fn)
    figure(fn);
    # Evaluate polynomial at derivatives:
    colorwheel = [:red,:blue,:green,:orange,:purple,:yellow,:cyan,:gray,:black];
    num_colors = 9;
    num_tsteps = 100;

    # points:
    p_inds = find(prob.poly_prob.B_orders.==0)
    px = prob.poly_prob.B_x[p_inds];
    py = prob.poly_prob.B_y[p_inds];

    subplot(2,2,1);
    scatter(px,py);

    for seg=1:prob.num_segs
        times = linspace(prob.times[seg],prob.times[seg+1],num_tsteps);
        degree = prob.degrees[seg]
        cx = prob.x_coeffs[:,seg];
        cy = prob.y_coeffs[:,seg];

        x0 = zeros(num_tsteps);
        x1 = zeros(num_tsteps);
        x2 = zeros(num_tsteps);
        x3 = zeros(num_tsteps);
        y0 = zeros(num_tsteps);
        y1 = zeros(num_tsteps);
        y2 = zeros(num_tsteps);
        y3 = zeros(num_tsteps);
        k=0;
        for t in times
            k+=1
            for d=1:degree
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
        # Plot the segment:
        subplot(2,2,1); title("Position");
        plot(x0,y0, color=colorwheel[mod(seg,num_colors)+1]);

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


function solve_polyseg_problem(prob::Poly2Segs)

    # Form A matrix:
    A = zeros(prob.tot_degree,prob.tot_degree);
    for constr=1:prob.tot_degree
        # Identify the polynomial:
        order = prob.poly_prob.B_orders[constr];
        t_ind = prob.poly_prob.B_time_inds[constr];
        time  = prob.times[t_ind];

        seg_ind = t_ind-1 # constraints always apply to the segment that ends at them
        if(seg_ind == 0)  # exception for initial constraints:
            seg_ind = 1
        end
        deg = prob.degrees[seg_ind];
        C_mat = prob.C_mats[:,:,seg_ind];

        # Now we compute the coefficient:
        for n=order:deg-1
            coeff=1;
            for k=1:order
                coeff *= (n+1-k);
            end
            pwr = n-order;
            A[constr,:] += coeff*(time^pwr).*C_mat[n+1,:]
        end
    end
     
    # Now solve for coefficients:
    x_reduced = A\prob.poly_prob.B_x;
    y_reduced = A\prob.poly_prob.B_y;

    for seg=1:prob.num_segs
        prob.x_coeffs[:,seg] = deepcopy(prob.C_mats[:,:,seg]*x_reduced);
        prob.y_coeffs[:,seg] = deepcopy(prob.C_mats[:,:,seg]*y_reduced);
    end
end


function gradient_descent(prob::Poly2Segs,step_size)
    num_points = prob.num_segs+1;

    # Compute cost for each segment by forming the Q matrices:
    costs_init = zeros(prob.num_segs) # vector of costs:
    costs_after = zeros(prob.num_segs) # vector of perturbed costs
    J=0;
    for seg=1:num_points-1
        Q_end = form_Q(prob.poly_prob.q_coeff, prob.times[seg+1]-prob.times[seg]); # Cost up to end point
        
        # Cost before perturbing: 
        costs_init[seg] = ((prob.x_coeffs[:,seg]'*Q_end*prob.x_coeffs[:,seg] + prob.y_coeffs[:,seg]'*Q_end*prob.y_coeffs[:,seg]))[1]
        J += costs_init[seg]
    end
    J+= prob.poly_prob.kT*prob.times[end]


    # Compute initial ratios:
    t_ratio = zeros(prob.num_segs);
    for seg=1:prob.num_segs
        t_ratio[seg] = (time_vec[seg+1]-time_vec[seg])/time_vec[end]
    end
    pert_size = 0.1; # perturb by 10%
    grad = zeros(prob.num_segs);
    pert_ratios = deepcopy(t_ratio);


    # Compute gradient for each segment:
    if(false)
    for seg= 1:prob.num_segs
        err = 0;
        for dir=-1:2:-1 # Optional to consider both increase an decrease. Default is to check decrease but otherwise increase.

            # Form new time vector based on perturbation
            t_new = zeros(prob.num_segs+1);
            t_new[1] = prob.times[1];        
            
            # Total increase in time:
            t_increase = (prob.times[seg+1]-prob.times[seg])*(dir*pert_size);
            # Relative decrease in everything else:
            decrease = prob.times[end]/(prob.times[end]+t_increase);                            

            for seg2=1:prob.num_segs
                ratio = 1.0;
                if(seg2==seg)
                    ratio=t_ratio[seg2]*(1+dir*pert_size)
                else
                    ratio=t_ratio[seg2]*decrease
                end
                dt = (prob.times[seg2+1]-prob.times[seg2])*ratio*prob.times[end];
                t_new[seg2+1] = dt + t_new[seg2];
            end
            
            # Compute gradient component:
            pnew = form_2segs(prob.poly_prob,t_new)
            solve_polyseg_problem(pnew)
            Q_end = form_Q(prob.poly_prob.q_coeff, t_new[seg+1]-t_new[seg])
            # Cost after perturbing: 
            costs_after[seg] = ((pnew.x_coeffs[:,seg]'*Q_end*pnew.x_coeffs[:,seg] + pnew.y_coeffs[:,seg]'*Q_end*pnew.y_coeffs[:,seg]))[1]

            # Compute cost for kT:
            costs_after[seg] += prob.poly_prob.kT*( (t_new[seg+1]-t_new[seg]))
            
            # Gradient is increase in cost for having moved.
            # This is dumb, but: 
            grad[seg] = (costs_after[seg]-costs_init[seg])/(max(costs_after[seg],costs_init[seg]))
            grad[seg] = sign(grad[seg])*pert_size;
#            if(costs_after[seg]-costs_init[seg] > 0) # Moving this way hurts, suggest going the other way:
#                grad[seg] += sign(dir)*(-1);
#            else
#                grad[seg] += sign(dir)
#            end
        end
        # This is dumb - fix later. Should depend on the error values. 
        # IF grad is positive, this means that both measurements agree that moving in the positive direction 
        if(grad[seg] == 0) # Means suggestion is to move positive
#            pert_ratios[pt] = (1-step_size)*t_ratio[pt] + (step_size)*(t_ratio[pt]*(1-pert_size)) # Go the other way
            println("Zero error!")
        else
            pert_ratios[seg] = (1-step_size)*t_ratio[seg] + step_size*(t_ratio[seg]*(1-grad[seg])) # move in suggested direction
        end
    end
    end

    println("1-G:\n",1-grad)
    println(pert_ratios)
    # Now compute gradient with respect to total time:
    t_scale = 1.0;
    t_new = zeros(prob.num_segs+1);
    t_new[1] = prob.times[1];
    for pt=2:num_points
        t_new[pt] = (1+pert_size)*prob.times[pt];
    end

    pnew = form_2segs(prob.poly_prob,t_new)
    solve_polyseg_problem(pnew)
    new_cost = 0;
    for seg=1:prob.num_segs
        Q_end = form_Q(prob.poly_prob.q_coeff, t_new[seg+1]-t_new[seg])
        # Cost after perturbing: 
        new_cost += ((pnew.x_coeffs[:,seg]'*Q_end*pnew.x_coeffs[:,seg] + pnew.y_coeffs[:,seg]'*Q_end*pnew.y_coeffs[:,seg]))[1]
    end
    new_cost += prob.poly_prob.kT*t_new[end];

    err = step_size*(new_cost-J)/max(new_cost,J)
    
    t_scale= 1-err;
    println("NC: $new_cost, J: $J");
    println("Tscale = $t_scale")

    if(J > new_cost) # increase helped
        t_scale = 1+pert_size;
    else  #Increase hurt
        t_scale = 1-pert_size;
    end    
    
    ## Now we take the gradient step:
    tvec_new = zeros(prob.num_segs+1);
    tvec_new[1] = prob.times[1];
    for seg=1:prob.num_segs
        dt = t_scale * prob.times[end]*pert_ratios[seg]
        tvec_new[seg+1] = dt + tvec_new[seg];
    end
    println("Tvec: $tvec_new");
    return tvec_new, J
end


function form_Abig(points, times, init_der, fin_der, cont_order)
    # Points and segments:
    num_pts = size(points,1);
    num_segs = num_pts-1;

    # Degrees:
    init_degree = size(init_der,1)+1+cont_order;
    fin_degree  = size(fin_der,1)+1+cont_order;
    other_degree = cont_order*2;
    seg_degrees = [init_degree; other_degree*ones(num_segs-2); fin_degree];
    seg_free_degrees = [size(init_der,1); ones(num_segs-2); size(fin_der,1)];

    tot_degree = init_degree+fin_degree+(num_segs-2)*other_degree;
    tot_free_degree = sum(seg_free_degrees);

    bval_fixed = [];

    # Form A matrix:
    A_big = zeros(tot_degree,tot_degree);
    a_ind = 1;


    num_fixed = size(init_der,1)+size(fin_der,1)+num_pts;
    num_free  = tot_free_degree;

    indep_degree = num_fixed + (cont_order-1)*(num_segs-1)


    c_ind_free = 1+num_fixed;
    c_ind_fixed = 1;
    c_ind = 1;
    C_big = zeros(tot_degree, indep_degree);

    for seg=1:num_segs
        curr_degree = round(Int64,seg_degrees[seg]);
        # arrange orders, degree into container to compute A_seg:
        orders_fixed = [];
        orders_free = [];
        times_fixed = [];
        times_free  = [];
        c_ind_init = c_ind;
        if(seg == 1)
            # Fixed derivatives:
            orders_fixed = [collect(0:size(init_der,1)); 0];
            nF = size(orders_fixed,1);
            C_big[c_ind:(c_ind+nF-1), c_ind_fixed:c_ind_fixed+nF-1] = eye(nF); # Fixed derivatives don't need to move
            c_ind_fixed += nF-1; # -1 b/c last point is redundant
            c_ind += nF;
            #         initial point, derivatives, end point
            bval_fixed  = [points[1]; init_der; points[2]];
#            times_fixed = [times[1]*ones(size(init_der,1)+1); times[2]];
            times_fixed = [zeros(size(init_der,1)+1); times[2]-times[1]];
            
            # Free Derivatives
            orders_free  = collect(1:cont_order-1); 
            nF = size(orders_free,1);
            C_big[(c_ind):(c_ind+nF-1), c_ind_free:c_ind_free+nF-1] = eye(nF); # Free derivatives given offset
            c_ind_free += nF - (cont_order-1)
            c_ind += nF
#            times_free  = times[2]*ones(cont_order-1);
            times_free = (times[2]-times[1])*ones(cont_order-1);

        elseif(seg == num_segs)
            orders_fixed = [0;collect(0:size(fin_der,1))];
            nF = size(orders_fixed,1);
            C_big[c_ind:(c_ind+nF-1), c_ind_fixed:c_ind_fixed+nF-1] = eye(nF); # Fixed derivatives don't need to move
            c_ind_fixed += nF;
            c_ind += nF;
            bval_fixed = [bval_fixed; points[end]; fin_der];
#            times_fixed = [times[end-1]; times[end]*ones(size(fin_der,1)+1)];
            times_fixed = [ 0; (times[end]-times[end-1])*ones(size(fin_der,1)+1)];

            orders_free =  collect(1:cont_order-1);
            nF = size(orders_free,1);
            C_big[(c_ind):(c_ind+nF-1), c_ind_free:c_ind_free+nF-1] = eye(nF); # Free derivatives given offset
            c_ind_free += nF - (cont_order-1)
            c_ind += nF
#            times_free  = times[end-1]*ones(cont_order-1);
            times_free = zeros(cont_order-1);
        else
            # Fixed derivatives
            orders_fixed = [0;0];
            bval_fixed  = [bval_fixed; points[seg+1]];
#            times_fixed = [times[seg-1]; times[seg]];
            times_fixed = [0; times[seg+1]-times[seg]];
            C_big[c_ind,c_ind_fixed] = 1; # Constraint on this segment start
            C_big[c_ind+1,c_ind_fixed+1] = 1; # Constraint on this segment end
            c_ind+=2; c_ind_fixed+=1;

            # Free derivatives:
            orders_free  = [collect(1:cont_order-1); collect(1:cont_order-1)];
#            times_free  = [times[seg-1]*ones(cont_order-1); times[seg]*ones(cont_order-1)];
            times_free = [zeros(cont_order-1); (times[seg+1]-times[seg])*ones(cont_order-1)];
            nF = size(orders_free,1);
            C_big[(c_ind):(c_ind+nF-1), c_ind_free:c_ind_free+nF-1] = eye(nF); # Free derivatives given offset
            c_ind_free += nF - (cont_order-1)
            c_ind += nF;
        end 

#        println("Segment $seg:");
#        println("Orders\t Times \t Values");
#        Cslice = C_big[c_ind_init:c_ind-1, 1:size(bval_fixed,1)]*bval_fixed;
#        for s = 1:size(orders_fixed,1)
#            println(orders_fixed[s], "\t", times_fixed[s], "\t",Cslice[s]); 
#        end


        A_seg = zeros(curr_degree,curr_degree);
        ords = [orders_fixed;orders_free];
        ts = [times_fixed;times_free];
        num_constr = size(ords,1);

        for k=1:num_constr
            A_seg[k,:] = constr_order(ords[k], ts[k],curr_degree);
        end

        A_big[a_ind:a_ind+curr_degree-1, a_ind:a_ind+curr_degree-1] = A_seg;
        a_ind += curr_degree;        
    end


    # Should also return C and B vectors

    return A_big,C_big,bval_fixed;
end

function form_Qbig(times, init_derivatives, final_derivatives, cont_order)
    if(size(init_derivatives,1) != cont_order-1)
        println("Warning - initial derivatives not the same length as continuity order")
    end
    if(size(final_derivatives,1) != cont_order-1)
        println("Warning - final derivatives not the same length as continuity order")
    end

    # Points and segments:
    num_pts = size(times,1);
    num_segs = num_pts-1;

    # Degrees:
    init_degree = size(init_derivatives,1)+cont_order+1;
    fin_degree  = size(final_derivatives,1)+cont_order+1;
    other_degree = cont_order*2;
    seg_degrees = [init_degree; other_degree*ones(num_segs-2); fin_degree];
    seg_free_degrees = [size(init_derivatives,1); ones(num_segs-2); size(final_derivatives,1)];

    tot_degree = init_degree+fin_degree+(num_segs-2)*other_degree;
    tot_free_degree = sum(seg_free_degrees);

    # Q matrix:
    Q = zeros(tot_degree,tot_degree);
    q_ind = 1;
    for seg=1:num_segs
        curr_degree = round(Int64, seg_degrees[seg]);
        q_coeffs = zeros(curr_degree); q_coeffs[4] = 1; q_coeffs[2] = 0.5; q_coeffs[3] = 0.5
        t = times[seg];
        if(seg!=1)
            t-=times[seg-1];
        end
        Q_i = form_Q(q_coeffs, t);
        Q[q_ind:q_ind+curr_degree-1, q_ind:q_ind+curr_degree-1] = Q_i;
        q_ind += curr_degree
    end
    return Q;    
end

