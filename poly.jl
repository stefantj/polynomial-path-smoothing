# For getting a handle on polygons and whatnot.
using PyPlot
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



num_points = 15
xvals = linspace(0,4*pi,num_points);
yvals = rand(num_points)*10
points = [xvals yvals zeros(num_points,2)];
figure(42); clf(); scatter(xvals, yvals);

vinit = [1,0,0,0]; vfinal =[1,0,0,0];
ainit = [0,0,0,0]; afinal =[0,0,0,0];

times = linspace(0,1,num_points);


# Unconstrained, single polynomial fitting:

degree = num_points + 4;
println("Running as a single $degree problem.");
# Coefficients
coeff_ux = zeros(degree); # these are our coefficients
coeff_uy = zeros(degree);

# Derivative constraint values
D_ux = [points[:,1]; vinit[1]; ainit[1]; vfinal[1]; afinal[1]];
D_uy = [points[:,2]; vinit[2]; ainit[2]; vfinal[2]; afinal[2]];


# Form constraint matrix A:

constr_orders = [zeros(num_points); 1;2;1;2];
constr_times  = [times; times[1]; times[1]; times[end]; times[end]];

# For this case, A_mat is full rank.
A_mat = zeros(degree,degree);
for constr=1:degree
    order = constr_orders[constr];
    time  = constr_times[constr];

    for n = order:degree-1
        coeff = 1;
        for k=1:order
            coeff *= (n+1-k);
        end
        pwr = n-order;
        A_mat[constr,n+1] = coeff * time^pwr;
    end
end

# Invert A matrix:
A_inv = inv(A_mat);
println("Condition number of single A: ", cond(A_inv));

# Back out coefficients:
coeff_ux = A_inv*D_ux;
coeff_uy = A_inv*D_uy;


# Plot
plot_poly(coeff_ux, coeff_uy, points, times[1], times[end], 42);


# Now Split into two segments, with three orders of continuity.
# Total degree is still 'degree', but each polynomial has a strange interaction - 
num_segs = 3;
initial_cont = 2;
final_cont = 2;
cont_order = 5; # order of continuity to enforce
# shuffle points around to balance order of segments
num_pointsi = round(Int64, (num_points+initial_cont+final_cont + (num_segs-1)*cont_order)/num_segs);
num_points2 = num_pointsi-cont_order;
num_points3 = num_pointsi-2-cont_order;
num_points1 = num_points - num_points2-num_points3;
println("There are 3 segments for $num_points: $num_points1+$num_points2+$num_points3");

indep_degree1 = 2+num_points1; # This is the number of constraints 1 will be responsible for
indep_degree2 = num_points2; # 
indep_degree3 = 2+num_points3; # 2 is for final points
degree = indep_degree1+indep_degree2+indep_degree3;

tot_degree1 = indep_degree1;
tot_degree2 = indep_degree2+cont_order
tot_degree3 = indep_degree3+cont_order

# Easy mapping
coeffs_map1 = [eye(indep_degree1) zeros(indep_degree1,degree-indep_degree1)];
coeffs_map2 = [zeros(cont_order, degree); [zeros(indep_degree2,indep_degree1) eye(indep_degree2) zeros(indep_degree2, indep_degree3)]];
coeffs_map3 = [zeros(cont_order, degree); [zeros(indep_degree3, indep_degree1+indep_degree2) eye(indep_degree3)]];

# Solve system backwards:
for k=cont_order:-1:1
    coeffs_map2[k,:] = coeffs_map1[k,:];
    max_deg = max(tot_degree1-1,tot_degree2-1)
    for i=k:max_deg
        coeff = times[num_points1]^(i-(k-1));
        for m=1:(k-1)
            coeff*=(i-(m-1))/(m);
        end
        if(i < tot_degree1)
            coeffs_map2[k,:] += coeff*coeffs_map1[i+1,:];
        end
        if(i < tot_degree2)
            coeffs_map2[k,:] -= coeff*coeffs_map2[i+1,:];
        end
    end
end

# Solve system backwards:
for k=cont_order:-1:1
    coeffs_map3[k,:] = coeffs_map2[k,:];
    max_deg = max(tot_degree3-1,tot_degree2-1)
    for i=k:max_deg
        coeff = times[num_points1+num_points2]^(i-(k-1));
        for m=1:(k-1)
            coeff*=(i-(m-1))/(m);
        end
        if(i < tot_degree2)
            coeffs_map3[k,:] += coeff*coeffs_map2[i+1,:];
        end
        if(i < tot_degree3)
            coeffs_map3[k,:] -= coeff*coeffs_map3[i+1,:];
        end
    end
end

# Derivative constraint values
D_ux = [points[:,1]; vinit[1]; ainit[1]; vfinal[1]; afinal[1]];
D_uy = [points[:,2]; vinit[2]; ainit[2]; vfinal[2]; afinal[2]];

# Now compute A matrix:

constr_orders = [zeros(num_points); 1;2;1;2];
constr_times  = [times; times[1]; times[1]; times[end]; times[end]];

# For this case, A_mat is full rank.
A_mat = zeros(degree,degree);
for constr=1:degree
    order = constr_orders[constr];
    time  = constr_times[constr];
    if(time <= times[num_points1]) # This constraint belongs to the first polynomial
        deg = tot_degree1;
        coeff_mat = coeffs_map1;
    elseif(time <= times[num_points1+num_points2]) # This constraint belongs to the second polynomial:
        deg = tot_degree2;
        coeff_mat = coeffs_map2;
    else # third.
        deg = tot_degree3;
        coeff_mat = coeffs_map3;
    end
    for n = order:deg-1
        coeff = 1;
        for k=1:order
            coeff *= (n+1-k);
        end
        pwr = n-order;
        A_mat[constr,:] += coeff*(time^pwr).*coeff_mat[n+1,:];
    end
end
# Invert A matrix:
A_inv = inv(A_mat);
println("Condition number of joint A: ", cond(A_mat)," ",cond(A_inv));

Ax = A_mat\D_ux;
Ay = A_mat\D_uy;
Ax = A_inv*D_ux;
Ay = A_inv*D_uy;

println("Errors: ", maximum(A_inv*A_mat-eye(degree)))

# Back out coefficients:
coeff_ux1 = coeffs_map1*Ax;
coeff_uy1 = coeffs_map1*Ay;
coeff_ux2 = coeffs_map2*Ax;
coeff_uy2 = coeffs_map2*Ay;
coeff_ux3 = coeffs_map3*Ax;
coeff_uy3 = coeffs_map3*Ay;

# Plot
plot_poly(coeff_ux1, coeff_uy1, points, times[1], times[num_points1], 43);
plot_poly(coeff_ux2, coeff_uy2, points, times[num_points1],times[num_points1+num_points2], 43);
plot_poly(coeff_ux3, coeff_uy3, points, times[num_points1+num_points2],times[end], 43);









