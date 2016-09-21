# Path smoothing based on Richter's work
#
using PyPlot
function plot_poly(cx,cy,pp, T, fn)
    figure(fn);
    radius = 50;
    degree = length(cx);
    times = linspace(0,T,100)
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
    
#    subplot(2,2,1); title("Position");
    plot(x0,y0, color=:green);
#    axis([-5,2*radius+5,-5,2*radius+5])
    axis([-5,65,-5,35])
    scatter(pp[:,1], pp[:,2]);
figure(5);
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


# Takes a vector of points and initial derivatives in
# spits a vector of differentially flat variables out (xyz \phi)
function smooth_path(point_path, v_init, a_init, v_fin, a_fin)

    # smooth paths in short segments if too long:
    max_pts = 3
    num_points = size(point_path,1)
    println("Planning path for $num_points points");

    v_i = v_init;
    a_i = a_init;
    v_f = nothing;
    a_f = nothing;

    t_segs = [];
    num_segs = 0;

    for curr_pt = 1:max_pts-1:num_points
        num_segs += 1;
        println("Segment $num_segs, starting at $curr_pt");
        if(curr_pt+max_pts-1 >= num_points)
            println("Last segment!, going to $num_points");
            v_f = v_fin;
            a_f = a_fin;
            pp = point_path[curr_pt:num_points,:];
        else
            println("going to ", curr_pt+max_pts-1)
            pp = point_path[curr_pt:(curr_pt+max_pts-1),:];
        end
        # Compute segment:
        xpoly,ypoly,zpoly,ppoly,T = smooth_short_path(pp, v_i, a_i, v_f, a_f);
        plot_poly(xpoly,ypoly,point_path,T,6)

        # Extract velocity, acceleration information
        v_i = zeros(4); a_i = zeros(4);
        for k = 1:length(xpoly)
            if(k > 1) 
                v_i[1] += (k-1)*(xpoly[k]*T^(k-2))
                v_i[2] += (k-1)*(ypoly[k]*T^(k-2))
                v_i[3] += (k-1)*(zpoly[k]*T^(k-2))
                v_i[4] += (k-1)*(ppoly[k]*T^(k-2))
            end
            if(k > 2)
                a_i[1] += (k-1)*(k-2)*(xpoly[k]*T^(k-3))
                a_i[2] += (k-1)*(k-2)*(ypoly[k]*T^(k-3))
                a_i[3] += (k-1)*(k-2)*(zpoly[k]*T^(k-3))
                a_i[4] += (k-1)*(k-2)*(ppoly[k]*T^(k-3))
            end
        end
        if(isnan(v_i[1]))
            v_i = zeros(4)
            a_i = zeros(4)
        end

        println("*** Initial velocity: ", norm(v_i), " Initial acc: ", norm(a_i))

        # Store coefficients for later
    end
end

function smooth_short_path(point_path, v_init, a_init, v_fin,a_fin)
    # Variables that describe the behavior of the smoother:
    num_points = size(point_path,1);   # Number of points in the path
    num_pts = num_points;
    num_segs = num_points-1;           # Number of segments in the path
    dof = 3;                           # Defines extra waypoints - weird things happen when this is not zero
    degree=0;
    if(v_fin == nothing)
        degree = (num_points+2) + dof
    else
        degree = (num_points+4) + dof      # Defines flexibility of the path - give it 
    end


    Q_coeffs = zeros(degree);
#    Q_coeffs[5] += .01 # Emphasize minimizing jerk
    Q_coeffs[4] += 1 # Emphasize minimizing jerk
#    Q_coeffs[3] += .01
#    Q_coeffs[2] += .01
#    Q_coeffs[1] += .01

    num_search_steps = 50;

    # perturbation
    t_diff = zeros(num_segs);
    J_curr= zeros(num_segs);
    J_last= zeros(num_segs);
    cost_trace = zeros(num_search_steps);

    best_J = Inf;
    best_x = zeros(degree);
    best_y = zeros(degree);
    best_z = zeros(degree);
    best_p = zeros(degree);
    best_T = zeros(num_points);

    

    # One polynomial for each flat variable
    coeff_x = zeros(degree);
    coeff_y = zeros(degree);
    coeff_z = zeros(degree);
    coeff_p = zeros(degree);

    # each polynomial is defined by degree derivatives
    D_x = zeros(degree);
    D_y = zeros(degree);
    D_z = zeros(degree);
    D_p = zeros(degree);

    for i = 1:num_segs
        D_x[i]=point_path[i+1,1]
        D_y[i]=point_path[i+1,2]
        D_z[i]=point_path[i+1,3]
        D_p[i]=point_path[i+1,4]
    end
    deg_ind = num_pts;
    if(v_fin==nothing)
        deg_ind -=2;
    else
        D_x[deg_ind] = v_fin[1]
        D_x[deg_ind+1] = a_fin[1]
    end
    D_x[deg_ind+2] = point_path[1,1];
    D_x[deg_ind+3] = v_init[1]
    D_x[deg_ind+4] = a_init[1]

    if(v_fin!=nothing)
        D_y[deg_ind] = v_fin[2]
        D_y[deg_ind+1] = a_fin[2]
    end
    D_y[deg_ind+2] = point_path[1,2];
    D_y[deg_ind+3] = v_init[2]
    D_y[deg_ind+4] = a_init[2]

    if(v_fin!=nothing)
        D_z[deg_ind] = v_fin[3]
        D_z[deg_ind+1] = a_fin[3]
    end
    D_z[deg_ind+2] = point_path[1,3];
    D_z[deg_ind+3] = v_init[3]
    D_z[deg_ind+4] = a_init[3]

    if(v_fin!=nothing)
        D_p[deg_ind] = v_fin[4]
        D_p[deg_ind+1] = a_fin[4]
    end
    D_p[deg_ind+2] = point_path[1,4];
    D_p[deg_ind+3] = v_init[4]
    D_p[deg_ind+4] = a_init[4]


    # This is where an inner loop could/should go to optimize T_split

    # Initial hueristic
    T_split = zeros(num_points);
    for k=1:num_segs
        T_split[k+1] = norm(point_path[k,:]-point_path[k+1,:])/norm(v_init)*(.5^k) + T_split[k]
    end


    step_size = .9

    T_min = 0.001; 
    T_max = T_split[end];

    for step=1:num_search_steps
    println("$step");

    tic();


    T_budget = T_split[end]

    if(T_budget > T_max)
        println(T_split)
#        println("Massive error!");
#        return;
    end


    println("T_split = $T_split");
    println("t_diff = $t_diff");
    # Form the constraints matrix.
    Q_mat = form_Q(Q_coeffs, T_budget);


    # Could change the order so A is upper triangular? might be a good idea.
    # Add constraints in the order they were put in:
    A_mat = zeros(degree,degree)

    constraint_times = [T_split[2:end]; T_split[end]; T_split[end]; zeros(3); zeros(dof)]
    constraint_order = [round(Int64, zeros(num_segs)); 1; 2;  round(Int64, (0:(dof+2)))]
    if(v_fin==nothing)
        # remove the num_pts,num_pts+1st elements
        constraint_times = [constraint_times[1:num_segs]; constraint_times[num_pts+2:end]];
        constraint_order = [constraint_order[1:num_segs]; constraint_order[num_pts+2:end]];
    end
    
    for constr=1:degree
    
        if(constraint_order[constr] >= degree) # Taking a derivative this large returns 0 always.
            continue;   # Zeros
        end
        order = constraint_order[constr];
        time  = constraint_times[constr];

        coeff = 1;
        # n is the coefficient
        for n = order:degree-1
            # zero if n < order
            coeff = 1;
            for k = 1:order
                coeff = coeff*(n+1-k);
            end 
        
            pwr = n - order;
            if(pwr <= 0)
                A_mat[constr,n+1] = coeff;
            else
                A_mat[constr,n+1] = coeff*(float(time)^(pwr)); 
            end
        end
    end

        
    t_init = toq();
    tic();

    # Now compute the R matrix.
    # Use block inverse identity:
    # M = [A B
    #      C 0]
    #
    # Minv = [0     Cinv
    #         Binv -Binv*A*Cinv] 


    
    n_A = num_segs+2;
    if(v_fin == nothing)
        n_A = num_segs
    end

    # This is down to microoptimization...
    # A can be whatever
    A_A = A_mat[1:n_A,1:(degree-n_A)]; 
    # B needs to be square
    A_B = A_mat[1:n_A, degree-n_A+1:degree];
    A_C = A_mat[n_A+1:degree, 1:(degree-n_A)];
    A_D = A_mat[n_A+1:degree, degree-n_A+1:degree];


    A_C_inv = diagm(1./diag(A_C))
    A_B_inv = eye(size(A_B,1));
    householder_Ldiv!(A_B, A_B_inv)
    A_inv = [A_D A_C_inv; A_B_inv -A_B_inv*A_A*A_C_inv];
    t_ainv = toq();

    tic();
    if(dof > 0)
        R = A_inv'*Q_mat*A_inv;    
        t_r1=toq();

        tic();

        # Compute gradient and set to zero
        mod = 4;
        if(v_fin == nothing)
            mod=2;
        end
        R_pp = R[num_points+mod+1:degree, num_points+mod+1:degree];
        R_fp = R[1:num_points+mod, num_points+mod+1:degree];

        opt_mat = -R_fp'; 
        householder_Ldiv!(R_pp, opt_mat);

        t_r_comp = toq();
        tic();
    # Compute values for the free derivatives
        D_x_opt = opt_mat*D_x[1:num_pts+mod];
        D_y_opt = opt_mat*D_y[1:num_pts+mod];
        D_z_opt = opt_mat*D_z[1:num_pts+mod];
        D_p_opt = opt_mat*D_p[1:num_pts+mod];
    t_precomp = toq();


    tic();
    # Reconstruct D vector:
        for i=1:(degree-num_pts-mod)
            D_x[num_pts+mod+i] = D_x_opt[i];
            D_y[num_pts+mod+i] = D_y_opt[i];
            D_z[num_pts+mod+i] = D_z_opt[i];
            D_p[num_pts+mod+i] = D_p_opt[i];
        end
    end

    # Now back out the polynomials
    coeff_x = A_inv*D_x;
    coeff_y = A_inv*D_y;
    coeff_z = A_inv*D_z;
    coeff_p = A_inv*D_p;
    t_finish = toq();
    t_qr = 0;
#    t_total = t_init+t_ainv+t_r1+t_r_comp+t_precomp+t_finish;
#    println("total compute time is $t_total. Init (", t_init/t_total, "), Ainv (", t_ainv/t_total, ") R (",t_r1/t_total,", ", t_r_comp/t_total ,") Pre (", t_precomp/t_total, ") Finish (", t_finish/t_total,")");


    # Compute gradient and move toward steepest descent
    # The optimization variable here is the _ratio_ between time per segments (minus the last one).
    for k = 1:num_segs
        J_last[k] = J_curr[k];
        J_curr[k] = 0;
        t = T_split[k+1];
        Q_mat = form_Q(Q_coeffs, t)

        J_curr[k] = (coeff_x'*Q_mat*coeff_x + coeff_y'*Q_mat*coeff_y + coeff_z'*Q_mat*coeff_z + coeff_p'*Q_mat*coeff_p)[1]
    end

    J_tot = sum(J_curr);
    cost_trace[step] = J_tot;
    J_last_tot = sum(J_last);
    # t_diff[1:num_segs-1] represents the change in ratio: There are num_segs-1 ratios
    # t_diff[num_segs] represents the change in total time
    
    println("Total cost: ", J_tot);

    if(J_tot < best_J)
        best_J = J_tot
        for k=1:degree
            best_x[k] = coeff_x[k];
            best_y[k] = coeff_y[k];
            best_z[k] = coeff_z[k];
            best_p[k] = coeff_p[k];
        end
        for k=1:num_points
            best_T[k] = T_split[k];
        end
    end


    if(step < num_search_steps)
        grad = zeros(num_segs)
        if(step != 1)
            for k = 1:num_segs-1
                grad[k] = (J_curr[k]-J_last[k]);#/log(t_diff[k]);
                if(grad[k] > 0)
                    if(t_diff[k] > 1)
                        t_diff[k] = .9
                    else
                        t_diff[k] = 1.1
                    end
                elseif(grad[k] < 0)
                    # don't change t_diff

                else
                    t_diff[k] = 1.0;
                end
    #            t_diff[k] = (grad[k])*step_size^step;
            end
            grad[num_segs] = (J_tot - J_last_tot)/t_diff[num_segs]
            t_diff[num_segs] = -grad[num_segs]*step_size^step

            println("Gradient: $grad");
        else
            # t_diff is arbitrary
            t_diff = rand(num_segs) + .5 # Multiply all variables by 50-150%
        end

        T_old = zeros(T_split)
        for k= 1:length(T_split)
            T_old[k] = T_split[k]
        end

        ratios = zeros(num_segs);
        for k=1:num_segs
            ratios[k] = (T_split[k+1]-T_split[k])/T_split[end];
        end
        println("Ratios $ratios");
        
        # Max time multiplies:
        if(T_split[end]+t_diff[end] < 0)
            t_diff[end] = -.1*T_split[end];
        end
        if(T_split[end]+t_diff[end] > T_max)
            t_diff[end] = T_max -T_split[end];
        end
        T_split[end] = T_split[end]+t_diff[end];
        norm_const = 0;
        for k = 1:num_segs-1
            norm_const += ratios[k]*t_diff[k]
            T_split[k+1] = T_split[k] + (ratios[k]*t_diff[k])*T_split[end]
        end
#        if(T_split[end] > T_max)
#            T_split./norm_const;
#        end
#        if(T_split[end] < T_min)
#            T_split .* T_min/T_split[end];
#        end
    end


    # Now if T_split is out of bounds, put back

        

#    plot_poly(coeff_x,coeff_y,point_path,T_split[end], 3)

    end #outer optimization loop

#    figure(4); plot(log(cost_trace))

#    plot_poly(best_x,best_y,point_path, best_T[end], 3)
    


    return best_x, best_y, best_z, best_p, best_T[end]
end

# Own implementation of householder algorithm
# Returns the result of (A^-1) M
# 
# A^-1M = R^-1 Q^T M
function householder_Ldiv!(A, M)
    m,n = size(A);
    R=A
    Res = M;
    m_m = size(M,2)
    println("A = "); println(A);
    println("B = "); println(M);

    for k=1:n
        y = R[k:m,k];
        v_k = y
        v_k[1] += sign(y[1])*norm(y)
        v_k /= norm(v_k)

        println("V_$k = $v_k");
        R[k:m, k:n] = R[k:m, k:n] - 2*v_k*(v_k'*R[k:m, k:n])

        for col = 1:m_m
            coeff = (v_k'*vec(Res[k:m,col]))[1];
            Res[k:m, col] -= 2*coeff*v_k;
        end
    end

    println("R = "); println(R);
    println("Q'B = ");println(Res);


    # Res now holds Q'M
    # Now we invert R and multiply. Result copied into Res
    backsubstitute!(R, Res);

    println("R\\B = "); println(Res);
end

# Own implementation of householder algorithm
# Returns the result of (A^-1) M
# 
# A^-1M = R^-1 Q^T M
function householder_Ldiv3!(A, M1, M2, M3)
    m,n = size(A);

    if(n > m)
        println("Error: A is the wrong shape (must be square or tall)");
        return
    end
    if(n < m)
        println("Warning: Inaccurate results for nonsquare A");
    end
    
    R = deepcopy(A);
    Res1 = M1;
    m_m1= size(M1,2)
    Res2 = M2;
    m_m2= size(M2,2)
    Res3 = M3;
    m_m3= size(M3,2)

    for k=1:n
        y = R[k:m,k];
        v_k = y
        v_k[1] += sign(y[1])*norm(y)
        v_k /= norm(v_k)

        R[k:m, k:n] = R[k:m, k:n] - 2*v_k*(v_k'*R[k:m, k:n])

        for col = 1:m_m1
            coeff = (v_k'*vec(Res1[k:m,col]))[1];
            Res1[k:m, col] -= 2*coeff*v_k;
        end
        for col = 1:m_m2
            coeff = (v_k'*vec(Res2[k:m,col]))[1];
            Res2[k:m, col] -= 2*coeff*v_k;
        end
        for col = 1:m_m3
            coeff = (v_k'*vec(Res3[k:m,col]))[1];
            Res3[k:m, col] -= 2*coeff*v_k;
        end
    end


    R = R[1:n, 1:n]

    # Res now holds Q'M
    # Now we invert R and multiply. Result copied into Res
    backsubstitute!(R, Res1);
    backsubstitute!(R, Res2);
    backsubstitute!(R, Res3);
end

# Backsubstitution to invert 
function backsubstitute!(M, Q)
    # M is an upper triangular matrix
    # use backsubstitution to invert efficiently.

    # For each column, run backsubstitution:
    num_rows,num_cols = size(Q);
    for col=1:num_cols
        for row=num_rows:-1:1
            rowsum=0;
            for k=num_rows:-1:row+1
                rowsum += Q[k,col]*M[row,k];
            end
            Q[row,col] = (Q[row,col]-rowsum)/M[row,row];
        end
    end
end

function form_Q(Q_coeffs, t)
    degree = length(Q_coeffs)
    Q_mat = zeros(degree,degree)

    for k = 0:degree-1
        if(Q_coeffs[k+1] == 0)
            continue;
        end
        c_k = Q_coeffs[k+1];
        for i=1:degree   
        # Form Q
        # Form Hessian. Assume minimizing 3rd derivatives (jerk, for vision)
            for l = 1:degree
                if(i >= min(k,l) && l >= k)
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


function test_householder()
    A = diagm(20:29);
    B = rand(10,10);
    Bcopy = zeros(10,10);
    for i=1:10
        for j=1:10
            B[i,j] = (i+19)+(j+19)*(j+19);
            Bcopy[i,j] = B[i,j]
        end
    end

    householder_Ldiv!(A,B);



end
