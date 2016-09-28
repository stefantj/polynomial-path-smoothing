using PyPlot # For plotting

type PolyProblem
    # Constraints
    B_x::Vector{Float64}        # Vector of constraint values for x dimension
    B_y::Vector{Float64}        # Vector of constraint values for y dimension
    B_z::Vector{Float64}        # Vector of constraint values for z dimension
    B_p::Vector{Float64}        # Vector of constraint values for phi dimension
    B_orders::Vector{Int64}     # Vector of orders of constraints
    B_time_inds::Vector{Int64}  # Vector of indices of time for constraints
    # Cost behavior
    q_coeff::Vector{Float64}    # Cost vector
    kT::Float64                 # Time cost
end

type PolyParams
    cont_order::Int64           # Order of continuity
end

type PolySol                  # For 1 segment per pair of points
    poly_prob::PolyProblem      # Original poly problem
    params::PolyParams
    num_segs::Int64             # Number of segments
    times::Vector{Float64}      # Times for each segment
    x_coeffs::Array{Float64,1}  # X coefficients
    y_coeffs::Array{Float64,1}  # Y coefficients
    z_coeffs::Array{Float64,1}  # X coefficients
    p_coeffs::Array{Float64,1}  # Y coefficients
end


# Forms a vanilla problem
function random_poly_prob(num_points::Int64, cont_order)
    # Random points
    xpts = sort(num_points*rand(num_points));
    ypts = num_points*rand(num_points);
    zpts = num_points*rand(num_points);
    ppts = num_points*rand(num_points);

    # Initial derivative constraints: 
    #          vel acc jerk
    bx_init = [1.0; zeros(cont_order-1)];
    by_init = [0.0; zeros(cont_order-1)];
    bz_init = zeros(cont_order);
    bp_init = zeros(cont_order);
    # Final derivative constraints:
    bx_final = [1.0; zeros(cont_order-1)];
    by_final = [0.0; zeros(cont_order-1)];
    bz_final = zeros(cont_order);
    bp_final = zeros(cont_order);

    #### Assemble constraint variables: #### 
    # Total B vectors:
    B_x = [xpts[1]; bx_init; xpts[2:end]; bx_final];
    B_y = [ypts[1]; by_init; ypts[2:end]; by_final];
    B_z = [zpts[1]; bz_init; zpts[2:end]; bz_final];
    B_p = [ppts[1]; bp_init; ppts[2:end]; bp_final];
    # Order of constraints, in order of B_x
    B_orders = [collect(0:cont_order-1); zeros(num_points-1); collect(1:cont_order-1)];
    # Time indices of constraints:
    B_time_inds = [ones(cont_order); collect(2:num_points); num_points*ones(cont_order-1)];

    degree = num_points+2*(cont_order-1);
    #### Define our cost metrics ####
    # Q coeffs: weight which derivatives you care about:
    q_coeff = zeros(degree);
    q_coeff[1] = 1.0;
    q_coeff[2] = 0.5;
    q_coeff[3] = 0.1;
    q_coeff[4] = 1.0;
    # Total time penalty
    kT = 50000; 

    return PolyProblem(B_x,B_y,B_z,B_p,round(Int64,B_orders),round(Int64,B_time_inds),q_coeff,kT)
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

function plot_poly(sol::PolySol, fn)
    # Evaluate polynomial at derivatives:
    colorwheel = [:red,:blue];
    num_colors = 2;
    num_tsteps = 100;


    # points:
    p_inds = find(sol.poly_prob.B_orders.==0)
    xpts = sol.poly_prob.B_x[p_inds];
    ypts = sol.poly_prob.B_y[p_inds];

    # plot polynomial:
    figure(fn);clf();
    offset=0;
    seg_ind = 1;
    for seg=1:sol.num_segs
        seg_deg = 0;
        if(seg==1)
            seg_deg = sol.params.cont_order + size(find(sol.poly_prob.B_time_inds.==1),1)
        elseif(seg==sol.num_segs)
            seg_deg = sol.params.cont_order + size(find(sol.poly_prob.B_time_inds.==sol.num_segs+1),1);
        else
            seg_deg = 2*sol.params.cont_order;
        end
        xc = sol.x_coeffs[seg_ind:seg_ind+seg_deg-1];
        yc = sol.y_coeffs[seg_ind:seg_ind+seg_deg-1];
        seg_ind+=seg_deg;
        
        t = linspace(0,sol.times[seg+1]-sol.times[seg],num_tsteps );
        xvals = zeros(num_tsteps);
        yvals = zeros(num_tsteps);
        xvels = zeros(num_tsteps);
        yvels = zeros(num_tsteps);
        xaccs = zeros(num_tsteps);
        yaccs = zeros(num_tsteps);
        xjerk = zeros(num_tsteps);
        yjerk = zeros(num_tsteps);

        for deg=1:seg_deg
            xvals += xc[deg]*t.^(deg-1);
            yvals += yc[deg]*t.^(deg-1);
        end

        for deg=2:seg_deg
            xvels += (deg-1)*xc[deg]*t.^(deg-2);
            yvels += (deg-1)*yc[deg]*t.^(deg-2);
        end

        for deg=3:seg_deg
            xaccs += (deg-1)*(deg-2)*xc[deg]*t.^(deg-3);
            yaccs += (deg-1)*(deg-2)*yc[deg]*t.^(deg-3);
        end

        for deg=4:seg_deg
            xjerk += (deg-1)*(deg-2)*(deg-3)*xc[deg]*t.^(deg-4);
            yjerk += (deg-1)*(deg-2)*(deg-3)*yc[deg]*t.^(deg-4);
        end
        figure(1);
        subplot(2,2,1); title("position");xlabel("X value"); ylabel("Y value");
        scatter(xpts,ypts); plot(xpts,ypts,color=:gray,linestyle=":");
        plot(xvals,yvals,color=colorwheel[mod(seg,num_colors)+1]);
        subplot(2,2,2); title("Velocity"); ylabel("Velocity"); xlabel("time");
        plot(t+sol.times[seg], xvels, color=:red);
        plot(t+sol.times[seg], yvels, color=:blue);
        legend(["X","Y"]);
        subplot(2,2,3); title("Acceleration"); ylabel("Acceleration"); xlabel("time");
        plot(t+sol.times[seg], xaccs, color=:red);
        plot(t+sol.times[seg], yaccs, color=:blue);
        legend(["X","Y"]);
        subplot(2,2,4); title("Jerk"); ylabel("Jerk"); xlabel("time");
        plot(t+sol.times[seg], xjerk, color=:red);
        plot(t+sol.times[seg], yjerk, color=:blue);
        legend(["X","Y"]);
    end
end


# Basic helper functions for poly plotting: 

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
function constr_order(order, time, degree)
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

# For multi-segment optimization
function form_Abig(prob::PolyProblem, param::PolyParams, times::Vector{Float64})

    init_time=0;
    prep_time=0; 
    addA_time=0;
tic();
    # Points and segments:
    num_pts = size(find(prob.B_orders.==0),1);
    num_segs = num_pts-1;

    cont_order = param.cont_order;

    # Degrees:
    num_init_constr = size(find(prob.B_time_inds.==1),1);
    num_fin_constr  = size(find(prob.B_time_inds.==num_pts),1);

    init_degree = num_init_constr+cont_order;
    fin_degree  = num_fin_constr+cont_order;

    seg_degrees = [init_degree; 2*cont_order*ones(num_segs-2); fin_degree];
    tot_degree = init_degree+fin_degree+(num_segs-2)*2*cont_order
    num_fixed = num_init_constr+num_fin_constr+num_pts-2;
    indep_degree = num_fixed + (cont_order-1)*(num_segs-1)

    # Form A matrix:
    A_big = zeros(tot_degree,tot_degree);
    A_big_inv = zeros(tot_degree,tot_degree);
    a_ind = 1;


    c_ind_free = 1+num_fixed;
    c_ind_fixed = 1;
    c_ind = 1;
    C_big = zeros(tot_degree, indep_degree);
init_time = toq();
    for seg=1:num_segs
        tic();
        curr_degree = round(Int64,seg_degrees[seg]);
        # arrange orders, degree into container to compute A_seg:
        orders_fixed = [];
        orders_free = [];
        times_fixed = [];
        times_free  = [];
        c_ind_init = c_ind;
        if(seg == 1)
            # Fixed derivatives:
            orders_fixed = [collect(0:num_init_constr-1); 0];
            nF = size(orders_fixed,1);
            C_big[c_ind:(c_ind+nF-1), c_ind_fixed:c_ind_fixed+nF-1] = eye(nF); # Fixed derivatives don't need to move
            c_ind_fixed += nF-1; # -1 b/c last point is redundant
            c_ind += nF;
            times_fixed = [zeros(num_init_constr); times[2]-times[1]];
            
            # Free Derivatives
            orders_free  = collect(1:cont_order-1); 
            nF = size(orders_free,1);
            C_big[(c_ind):(c_ind+nF-1), c_ind_free:c_ind_free+nF-1] = eye(nF); # Free derivatives given offset
            c_ind_free += nF - (cont_order-1)
            c_ind += nF
            times_free = (times[2]-times[1])*ones(cont_order-1);

        elseif(seg == num_segs)
            orders_fixed = [0;collect(0:num_fin_constr-1)];
            nF = size(orders_fixed,1);
            C_big[c_ind:(c_ind+nF-1), c_ind_fixed:c_ind_fixed+nF-1] = eye(nF); # Fixed derivatives don't need to move
            c_ind_fixed += nF;
            c_ind += nF;
            times_fixed = [ 0; (times[end]-times[end-1])*ones(num_fin_constr)];

            orders_free =  collect(1:cont_order-1);
            nF = size(orders_free,1);
            C_big[(c_ind):(c_ind+nF-1), c_ind_free:c_ind_free+nF-1] = eye(nF); # Free derivatives given offset
            c_ind_free += nF - (cont_order-1)
            c_ind += nF
            times_free = zeros(cont_order-1);
        else # Could precompute most of this
            # Fixed derivatives
            orders_fixed = [0;0];
            times_fixed = [0; times[seg+1]-times[seg]];
            C_big[c_ind,c_ind_fixed] = 1; # Constraint on this segment start
            C_big[c_ind+1,c_ind_fixed+1] = 1; # Constraint on this segment end
            c_ind+=2; c_ind_fixed+=1;

            # Free derivatives:
            orders_free  = [collect(1:cont_order-1); collect(1:cont_order-1)];
            times_free = [zeros(cont_order-1); (times[seg+1]-times[seg])*ones(cont_order-1)];
            nF = size(orders_free,1);
            C_big[(c_ind):(c_ind+nF-1), c_ind_free:c_ind_free+nF-1] = eye(nF); # Free derivatives given offset
            c_ind_free += nF - (cont_order-1)
            c_ind += nF;
        end 
        prep_time+=toq();
        tic();
        A_seg = zeros(curr_degree,curr_degree);
        ords = [orders_fixed;orders_free];
        ts = [times_fixed;times_free];
        num_constr = size(ords,1);

        for k=1:num_constr
            A_seg[k,:] = constr_order(ords[k], ts[k],curr_degree);
        end

        A_big[a_ind:a_ind+curr_degree-1, a_ind:a_ind+curr_degree-1] = A_seg;

    # Since A_big is block diagonal, the inverse is the blocks of Abig
    # The blocks have a special structure:
    # [A B
    #  C 0]
    # A is upper triangular, B is full, C is diagonal and 'D' is 0 
    # some massaging is needed to make this work, but this is pretty easy to do.
    # Then we can compute the inverse of the very small 5x5, which will speed things up a lot.
        A_big_inv[a_ind:a_ind+curr_degree-1, a_ind:a_ind+curr_degree-1] = inv(A_seg);
        a_ind += curr_degree;        
        addA_time+=toq();
    end

    tot_time = init_time+prep_time+addA_time;
#    println("A_init: ", init_time/tot_time, "\nA_prep: ", prep_time/tot_time, "\nA_add: ", addA_time/tot_time);



    return A_big,A_big_inv,C_big;
end

function form_Qbig(times,num_init, num_fin, cont_order)
    if(num_init != cont_order)
        println("Warning - initial derivatives not the same length as continuity order")
    end
    if(num_fin != cont_order)
        println("Warning - final derivatives not the same length as continuity order")
    end

    # Points and segments:
    num_pts = size(times,1);
    num_segs = num_pts-1;

    # Degrees:
    init_degree = num_init+cont_order;
    fin_degree  = num_fin+cont_order;
    other_degree = cont_order*2;
    seg_degrees = [init_degree; other_degree*ones(num_segs-2); fin_degree];
    seg_free_degrees = [num_init-1; ones(num_segs-2); num_fin-1];

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


#function gradient_descent(prob::Poly2Segs,step_size)
#    num_points = prob.num_segs+1;

    # Compute cost for each segment by forming the Q matrices:
#    costs_init = zeros(prob.num_segs) # vector of costs:
#    costs_after = zeros(prob.num_segs) # vector of perturbed costs
#    J=0;
#    for seg=1:num_points-1
#        Q_end = form_Q(prob.poly_prob.q_coeff, prob.times[seg+1]-prob.times[seg]); # Cost up to end point
        
        # Cost before perturbing: 
#        costs_init[seg] = ((prob.x_coeffs[:,seg]'*Q_end*prob.x_coeffs[:,seg] + prob.y_coeffs[:,seg]'*Q_end*prob.y_coeffs[:,seg]))[1]
#        J += costs_init[seg]
#    end
#    J+= prob.poly_prob.kT*prob.times[end]


    # Compute initial ratios:
#    t_ratio = zeros(prob.num_segs);
#    for seg=1:prob.num_segs
#        t_ratio[seg] = (time_vec[seg+1]-time_vec[seg])/time_vec[end]
#    end
#    pert_size = 0.1; # perturb by 10%
#    grad = zeros(prob.num_segs);
#    pert_ratios = deepcopy(t_ratio);
#
#
#    # Compute gradient for each segment:
#    if(false)
#    for seg= 1:prob.num_segs
#        err = 0;
#        for dir=-1:2:-1 # Optional to consider both increase an decrease. Default is to check decrease but otherwise increase.
#
#            # Form new time vector based on perturbation
#            t_new = zeros(prob.num_segs+1);
#            t_new[1] = prob.times[1];        
#            
#            # Total increase in time:
#            t_increase = (prob.times[seg+1]-prob.times[seg])*(dir*pert_size);
#            # Relative decrease in everything else:
#            decrease = prob.times[end]/(prob.times[end]+t_increase);                            
#
#            for seg2=1:prob.num_segs
#                ratio = 1.0;
#                if(seg2==seg)
#                    ratio=t_ratio[seg2]*(1+dir*pert_size)
#                else
#                    ratio=t_ratio[seg2]*decrease
#                end
#                dt = (prob.times[seg2+1]-prob.times[seg2])*ratio*prob.times[end];
#                t_new[seg2+1] = dt + t_new[seg2];
#            end
#            
#            # Compute gradient component:
#            pnew = form_2segs(prob.poly_prob,t_new)
#            solve_polyseg_problem(pnew)
#            Q_end = form_Q(prob.poly_prob.q_coeff, t_new[seg+1]-t_new[seg])
#            # Cost after perturbing: 
#            costs_after[seg] = ((pnew.x_coeffs[:,seg]'*Q_end*pnew.x_coeffs[:,seg] + pnew.y_coeffs[:,seg]'*Q_end*pnew.y_coeffs[:,seg]))[1]
#
#            # Compute cost for kT:
#            costs_after[seg] += prob.poly_prob.kT*( (t_new[seg+1]-t_new[seg]))
#            
#            # Gradient is increase in cost for having moved.
#            # This is dumb, but: 
#            grad[seg] = (costs_after[seg]-costs_init[seg])/(max(costs_after[seg],costs_init[seg]))
#            grad[seg] = sign(grad[seg])*pert_size;
##            if(costs_after[seg]-costs_init[seg] > 0) # Moving this way hurts, suggest going the other way:
##                grad[seg] += sign(dir)*(-1);
##            else
##                grad[seg] += sign(dir)
##            end
#        end
#        # This is dumb - fix later. Should depend on the error values. 
#        # IF grad is positive, this means that both measurements agree that moving in the positive direction 
#        if(grad[seg] == 0) # Means suggestion is to move positive
##            pert_ratios[pt] = (1-step_size)*t_ratio[pt] + (step_size)*(t_ratio[pt]*(1-pert_size)) # Go the other way
#            println("Zero error!")
#        else
#            pert_ratios[seg] = (1-step_size)*t_ratio[seg] + step_size*(t_ratio[seg]*(1-grad[seg])) # move in suggested direction
#        end
#    end
#    end

#    println("1-G:\n",1-grad)
#    println(pert_ratios)
#    # Now compute gradient with respect to total time:
#    t_scale = 1.0;
#    t_new = zeros(prob.num_segs+1);
#    t_new[1] = prob.times[1];
#    for pt=2:num_points
#        t_new[pt] = (1+pert_size)*prob.times[pt];
#    end

#    pnew = form_2segs(prob.poly_prob,t_new)
#    solve_polyseg_problem(pnew)
#    new_cost = 0;
#    for seg=1:prob.num_segs
#        Q_end = form_Q(prob.poly_prob.q_coeff, t_new[seg+1]-t_new[seg])
#        # Cost after perturbing: 
#        new_cost += ((pnew.x_coeffs[:,seg]'*Q_end*pnew.x_coeffs[:,seg] + pnew.y_coeffs[:,seg]'*Q_end*pnew.y_coeffs[:,seg]))[1]
#    end
#    new_cost += prob.poly_prob.kT*t_new[end];
#
#    err = step_size*(new_cost-J)/max(new_cost,J)
#    
#    t_scale= 1-err;
#    println("NC: $new_cost, J: $J");
#    println("Tscale = $t_scale")
#
#    if(J > new_cost) # increase helped
#        t_scale = 1+pert_size;
#    else  #Increase hurt
#        t_scale = 1-pert_size;
#    end    
#    
#    ## Now we take the gradient step:
#    tvec_new = zeros(prob.num_segs+1);
#    tvec_new[1] = prob.times[1];
#    for seg=1:prob.num_segs
#        dt = t_scale * prob.times[end]*pert_ratios[seg]
#        tvec_new[seg+1] = dt + tvec_new[seg];
#    end
#    println("Tvec: $tvec_new");
#    return tvec_new, J
#end


# Smooths the given polynomial problem. 
# Answer is returned as a piecewise-polynomial with specified degree of continuity
function poly_smoothing(prob::PolyProblem, param::PolyParams)
    # Extract information to form Abig:
    num_points = maximum(prob.B_time_inds);

    num_fixed = size(prob.B_x,1);

    num_init_constr = size(find(prob.B_time_inds.==1),1);
    num_fin_constr  = size(find(prob.B_time_inds.==num_points),1);

    times = float(collect(0:num_points-1));

## Here is where the gradient loop will start:
    # Form A matrix:

    tic();
    A,Ainv,C = form_Abig(prob, param, times);
    num_unique = size(C,2);
    t_Abig = toq();

#    figure(5); spy(A);

    tic();
#    Ainv = inv(A);
    AiC = Ainv*C; # This is about 10% time, and can be fixed by just selecting rows/columns
    t_Ainv = toq();

    tic();
    # Form Q matrix:
    Q = form_Qbig(times, num_init_constr, num_fin_constr, param.cont_order); 
    t_Q = toq();

    tic();
    # Compute free variables:
    R = AiC'*(Q*AiC);
# should be a householder multiply here
    opt_mat = - ( R[num_fixed+1:num_unique, num_fixed+1:num_unique])\R[1:num_fixed, num_fixed+1:num_unique]';

    bx = [prob.B_x; opt_mat*prob.B_x];
    by = [prob.B_y; opt_mat*prob.B_y];
    bz = [prob.B_z; opt_mat*prob.B_z];
    bp = [prob.B_p; opt_mat*prob.B_p];

    x_coeffs = AiC*bx;
    y_coeffs = AiC*by;
    z_coeffs = AiC*bz;
    p_coeffs = AiC*bp;
    t_opt = toq();

    t_tot = t_Abig+t_Ainv+t_Q+t_opt;
#    println("Timing breakdown:");
#    println("Abig: ", t_Abig, "\t", t_Abig/t_tot);
#    println("Ainv: ", t_Ainv, "\t", t_Ainv/t_tot);
#    println("Q:    ", t_Q, "\t", t_Q/t_tot);
#    println("opt:  ", t_opt, "\t", t_opt/t_tot); 

## This is where the gradient should be updated and times recomputed

## Here is where the gradient loop stops

    sol = PolySol(prob, param, num_points-1, times, x_coeffs, y_coeffs, z_coeffs, p_coeffs); 
    return sol
end


function time_multiseg()
    # Make random problem:
    t1s = zeros(48);
    t2s = zeros(48);
    tind = 0;
    numtries = 10
    for num_points = 3:50
println("$num_points points");
            tind+=1;
        for t=1:numtries
            cont_order = 5;
            tic();
            prob = random_poly_prob(num_points,cont_order);
            t1s[tind]+=(toq()/numtries);
            params = PolyParams(cont_order);
            tic();
            sol = poly_smoothing(prob, params);
            t2s[tind]+=(toq()/numtries)
        end
    end

    figure(2); plot(collect(3:50), t1s,color=:red);
    plot(collect(3:50), t2s,color=:blue);
    ylabel("Computation time");
    xlabel("Number of points");
    legend(["Problem formation","Problem solution"]);
end

function test_multiseg(num_points)
    cont_order = 5;
    prob = random_poly_prob(num_points,cont_order);
    params = PolyParams(cont_order);
    sol = poly_smoothing(prob, params)
    plot_poly(sol,1);
end
