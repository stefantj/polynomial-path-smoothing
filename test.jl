# Code for testing the connect_poly code.

using PyPlot


include("poly_helper.jl");

# The order of the constraints is defined by the size of the vectors
#function connect_points(init_config::Vector{Point}, final_config::Vector{Point}, Q_coeffs)## >> These are already implemented/ Stefan can update them << ##
    # Make a polyseg to test. Note that this is not a real polysegment, and inconsistent.
#    x_coeffs = [1,1,1,-0.4,0];
#    y_coeffs = [1,3,-.3,.01,0];
#    z_coeffs = [0,1,0,0,0];
#    p_coeffs = zeros(5);
#    t=3.0
#    q = 0.0;
#    cells = [1,2,3,4,5,104,203,301,402];
#    init_config = [Point(0.0,0.0,0.0,0.0)];
#    final_config = [Point(0.0,0.0,0.0,0.0)];
#    return poly_segment(x_coeffs,y_coeffs,z_coeffs,p_coeffs,t,q,cells,init_config,final_config);
#end

##### Occupancy grid helper functions #### 
# Let's work in a 10m x 10m environment with resolution .1m
function get_grid_extent()
    return 10;
end
function get_grid_resolution()
    return 0.1
end

# function which returns a boolean whether the cell at (xyz) is occupied (true = occupied)
function occupancy_check(x::Float64,y::Float64,z::Float64)
    id = occupancy_get_id(x,y,z);
    return occupancy_check(id);
end

# Function which returns a boolean whether cell cell_id is occupied (true = occupied)
function occupancy_check(cell_id::Int64)
    return false; # No obstacles for now...
end

# Function which returns the cell id of the point (xyz). Returns an integer.
function occupancy_get_id(x,y,z)
    width = get_grid_extent();
    res   = get_grid_resolution();
    n = round(Int64,ceil(width/res));
    return sub2ind((n,n,n),round(Int64,x/res)+1,round(Int64,y/res)+1,round(Int64,z/res)+1)
end

function test_on_circle()

    # put in points as a circle:
    times = collect(0:0.01:2pi)
    xs = 3*sin(times)+6;
    ys = cos(times)+5;
    zs = 2*sin(times)+3;
    tind = 0;
   cells = Int64[];
    for t in times
        tind+=1;
        push!(cells, occupancy_get_id(xs[tind],ys[tind],zs[tind]))
    end 

    N = round(Int64, ceil(get_grid_extent()/get_grid_resolution())); # Number of points in a dimension
    xy_check_grid = zeros(N,N); # 0 means not checked.
    xz_check_grid = zeros(N,N); # 0 means not checked.
    yz_check_grid = zeros(N,N); # 0 means not checked.


    for cell in cells
        # Convert to subindices:  
        subs = ind2sub((N,N,N), cell);
        # Indicate that we check:
        xy_check_grid[subs[1],subs[2]] += 1;
        xz_check_grid[subs[1],subs[3]] += 1;
        yz_check_grid[subs[2],subs[3]] += 1;
    end 

    num_tsteps=100;                                                         

    width = get_grid_extent();
    delta = get_grid_resolution()/2;
    figure(1,figsize=(6,6)); clf();    
    subplot(2,2,1);
    imshow(xy_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
    plot(xs+delta,width-ys-delta,color=:red);    
    title("XY plane");
    subplot(2,2,2);
    imshow(yz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
    plot(ys+delta,width-zs-delta,color=:red)
    title("YZ plane");
    subplot(2,2,3);
    imshow(xz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
    plot(xs+delta,width-zs-delta,color=:red)
    title("XZ plane");

end

# Now we test:
function check_poly_segment(poly::poly_segment)
    # Try 3D plot
    width = get_grid_extent();
    delta = get_grid_resolution()/2;
    num_tsteps=100;
    times = linspace(0,poly.t,num_tsteps);
    xs = zeros(num_tsteps)    
    ys = zeros(num_tsteps);
    zs = zeros(num_tsteps);

    order = size(poly.x_coeffs,1);
    for ord = 1:order
        xs += poly.x_coeffs[ord].*(times.^(ord-1))
        ys += poly.y_coeffs[ord].*(times.^(ord-1))
        zs += poly.z_coeffs[ord].*(times.^(ord-1))
    end
    figure(3); clf();
    plot3D(xs,ys,zs);
    xlabel("x");
    ylabel("y");
    zlabel("z");
    #xlim([0,width]);
    #ylim([0,width]);
    #zlim([0,width]);
    
    
    
    
    #############################3Check configurations##########################
    # Should add more checks such as confirming velocity
    #Check if all configurations are within error between given and actual
    precision = 0.00000000001;
    #Check if there is any mismatch in configurations
    println("All initial configs good to within precision? If not, a Fail will appear here.")
    for g = 1:3;
        #initial condition check
        if(!(abs(poly.init_config[g].x - evaluate_poly(poly.x_coeffs, g-1, 0) <= precision)))
            println("Fail in the initial x ", g-1, "th deriv: ", poly.init_config[g].x, ", ", evaluate_poly(poly.x_coeffs, g-1, 0))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.x_coeffs)
        end
        if(!(abs(poly.init_config[g].y - evaluate_poly(poly.y_coeffs, g-1, 0) <= precision)))
            println("Fail in the initial y ", g-1, "th deriv: ", poly.init_config[g].y, ", ", evaluate_poly(poly.y_coeffs, g-1, 0))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.y_coeffs)
        end
        if(!(abs(poly.init_config[g].z - evaluate_poly(poly.z_coeffs, g-1, 0) <= precision)))
            println("Fail in the initial z ", g-1, "th deriv: ", poly.init_config[g].z, ", ", evaluate_poly(poly.z_coeffs, g-1, 0))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.z_coeffs)
        end        
        if(!(abs(poly.init_config[g].p - evaluate_poly(poly.p_coeffs, g-1, 0) <= precision)))
            println("Fail in theinitial p ", g-1, "th deriv: ", poly.init_config[g].p, ", ", evaluate_poly(poly.p_coeffs, g-1, 0))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.p_coeffs)        
        end
        #final condition check
        if(!(abs(poly.final_config[g].x - evaluate_poly(poly.x_coeffs, g-1, poly.t) <= precision)))
            println("Fail in the final x ", g-1, "th deriv: ", poly.final_config[g].x, ", ", evaluate_poly(poly.x_coeffs, g-1, poly.t))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.x_coeffs)            
        end
        if(!(abs(poly.final_config[g].y - evaluate_poly(poly.y_coeffs, g-1, poly.t) <= precision)))
            println("Fail in the final y ", g-1, "th deriv: ", poly.final_config[g].y, ", ", evaluate_poly(poly.y_coeffs, g-1, poly.t))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.y_coeffs)            
        end
        if(!(abs(poly.final_config[g].z - evaluate_poly(poly.z_coeffs, g-1, poly.t) <= precision)))
            println("Fail in the final z ", g-1, "th deriv: ", poly.final_config[g].z, ", ", evaluate_poly(poly.z_coeffs, g-1, poly.t))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.z_coeffs)            
        end        
        if(!(abs(poly.final_config[g].p - evaluate_poly(poly.p_coeffs, g-1, poly.t) <= precision)))
            println("Fail in the final p ", g-1, "th deriv: ", poly.final_config[g].p, ", ", evaluate_poly(poly.p_coeffs, g-1, poly.t))
            println(poly.init_config)
            println(poly.final_config)
            println(poly.y_coeffs)            
        end
    end
            
    #println(abs(poly.init_config[1].x - evaluate_poly(poly.x_coeffs, 0, 0)) < precision && 
    #abs(poly.init_config[1].y-evaluate_poly(poly.y_coeffs, 0, 0)) < precision &&
    #abs(poly.init_config[1].z-evaluate_poly(poly.z_coeffs, 0, 0)) < precision &&
    #abs(poly.init_config[1].p-evaluate_poly(poly.p_coeffs, 0, 0)) < precision &&
    #abs(poly.init_config[2].x-evaluate_poly(poly.x_coeffs, 1, 0)) < precision && 
    #abs(poly.init_config[2].y-evaluate_poly(poly.y_coeffs, 1, 0)) < precision &&
    #abs(poly.init_config[2].z-evaluate_poly(poly.z_coeffs, 1, 0)) < precision &&
    #abs(poly.init_config[2].p-evaluate_poly(poly.p_coeffs, 1, 0)) < precision &&
    #abs(poly.init_config[3].x-evaluate_poly(poly.x_coeffs, 2, 0)) < precision && 
    #abs(poly.init_config[3].y-evaluate_poly(poly.y_coeffs, 2, 0)) < precision &&
    #abs(poly.init_config[3].z-evaluate_poly(poly.z_coeffs, 2, 0)) < precision &&
    #abs(poly.init_config[3].p-evaluate_poly(poly.p_coeffs, 2, 0)) < precision &&
    #abs(poly.final_config[1].x-evaluate_poly(poly.x_coeffs, 0, poly.t)) < precision && 
    #abs(poly.final_config[1].y-evaluate_poly(poly.y_coeffs, 0, poly.t)) < precision &&
    #abs(poly.final_config[1].z-evaluate_poly(poly.z_coeffs, 0, poly.t)) < precision &&
    #abs(poly.final_config[1].p-evaluate_poly(poly.p_coeffs, 0, poly.t)) < precision &&
    #abs(poly.final_config[2].x-evaluate_poly(poly.x_coeffs, 1, poly.t)) < precision && 
    #abs(poly.final_config[2].y-evaluate_poly(poly.y_coeffs, 1, poly.t)) < precision &&
    #abs(poly.final_config[2].z-evaluate_poly(poly.z_coeffs, 1, poly.t)) < precision &&
    #abs(poly.final_config[2].p-evaluate_poly(poly.p_coeffs, 1, poly.t)) < precision &&
    #abs(poly.final_config[3].x-evaluate_poly(poly.x_coeffs, 2, poly.t)) < precision && 
    #abs(poly.final_config[3].y-evaluate_poly(poly.y_coeffs, 2, poly.t)) < precision &&
    #abs(poly.final_config[3].z-evaluate_poly(poly.z_coeffs, 2, poly.t)) < precision &&
    #abs(poly.final_config[3].p-evaluate_poly(poly.p_coeffs, 2, poly.t)) < precision
    #)   
    
    #println("Initial")
    #println(abs(poly.init_config[1].x - evaluate_poly(poly.x_coeffs, 0, 0)) < precision, ", ", 
    #abs(poly.init_config[1].y-evaluate_poly(poly.y_coeffs, 0, 0)) < precision  , ", ",
    #abs(poly.init_config[1].z-evaluate_poly(poly.z_coeffs, 0, 0)) < precision, ", ",
    #abs(poly.init_config[1].p-evaluate_poly(poly.p_coeffs, 0, 0)) < precision)
    #println(abs(poly.init_config[2].x - evaluate_poly(poly.x_coeffs, 1, 0)) < precision, ", ", 
    #abs(poly.init_config[2].y-evaluate_poly(poly.y_coeffs, 1, 0)) < precision  , ", ",
    #abs(poly.init_config[2].z-evaluate_poly(poly.z_coeffs, 1, 0)) < precision, ", ",
    #abs(poly.init_config[2].p-evaluate_poly(poly.p_coeffs, 1, 0)) < precision)
    #println(abs(poly.init_config[3].x - evaluate_poly(poly.x_coeffs, 2, 0)) < precision, ", ", 
    #abs(poly.init_config[3].y-evaluate_poly(poly.y_coeffs, 2, 0)) < precision  , ", ",
    #abs(poly.init_config[3].z-evaluate_poly(poly.z_coeffs, 2, 0)) < precision, ", ",
    #abs(poly.init_config[3].p-evaluate_poly(poly.p_coeffs, 2, 0)) < precision)
    #println("Final")
    #println(abs(poly.final_config[1].x - evaluate_poly(poly.x_coeffs, 0, poly.t)) < precision, ", ", 
    #abs(poly.final_config[1].y-evaluate_poly(poly.y_coeffs, 0, poly.t)) < precision  , ", ",
    #abs(poly.final_config[1].z-evaluate_poly(poly.z_coeffs, 0, poly.t)) < precision, ", ",
    #abs(poly.final_config[1].p-evaluate_poly(poly.p_coeffs, 0, poly.t)) < precision)
    #println(abs(poly.final_config[2].x - evaluate_poly(poly.x_coeffs, 1, poly.t)) < precision, ", ", 
    #abs(poly.final_config[2].y-evaluate_poly(poly.y_coeffs, 1, poly.t)) < precision  , ", ",
    #abs(poly.final_config[2].z-evaluate_poly(poly.z_coeffs, 1, poly.t)) < precision, ", ",
    #abs(poly.final_config[2].p-evaluate_poly(poly.p_coeffs, 1, poly.t)) < precision)
    #println(abs(poly.final_config[3].x - evaluate_poly(poly.x_coeffs, 2, poly.t)) < precision, ", ", 
    #abs(poly.final_config[3].y-evaluate_poly(poly.y_coeffs, 2, poly.t)) < precision  , ", ",
    #abs(poly.final_config[3].z-evaluate_poly(poly.z_coeffs, 2, poly.t)) < precision, ", ",
    #abs(poly.final_config[3].p-evaluate_poly(poly.p_coeffs, 2, poly.t)) < precision)
    #println("The given initial points")
    #println(poly.init_config[1].x, ", ", poly.init_config[1].y, ", ", poly.init_config[1].z, ", ", poly.init_config[1].p)
    #Print the configurations and compare to what is calculated within error
    #println("The initial points of the polynomial")
    #println(evaluate_poly(poly.x_coeffs, 0, 0), ", ", evaluate_poly(poly.y_coeffs, 0, 0), ", ", evaluate_poly(poly.z_coeffs, 0, 0), ", ", evaluate_poly(poly.p_coeffs, 0, 0))
    #initial velocities
    #println("The given initial velocities")
    #println(poly.init_config[2].x, ", ", poly.init_config[2].y, ", ", poly.init_config[2].z, ", ", poly.init_config[2].p)
    #Print the configurations and compare to what is calculated within error
    #println("The initial velocities of the polynomial")
    #println(evaluate_poly(poly.x_coeffs, 1, 0), ", ", evaluate_poly(poly.y_coeffs, 1, 0), ", ", evaluate_poly(poly.z_coeffs, 1, 0), ", ", evaluate_poly(poly.p_coeffs, 1, 0))
    #initial accelerations
    #println("The given initial accelerations")
    #println(poly.init_config[3].x, ", ", poly.init_config[3].y, ", ", poly.init_config[3].z, ", ", poly.init_config[3].p)
    #Print the configurations and compare to what is calculated within error
    #println("The initial accelerations of the polynomial")
    #println(evaluate_poly(poly.x_coeffs, 2, 0), ", ", evaluate_poly(poly.y_coeffs, 2, 0), ", ", evaluate_poly(poly.z_coeffs, 2, 0), ", ", evaluate_poly(poly.p_coeffs, 2, 0))
    
    
    ###########Print the final configuarations and how they compare to what the polynomial gives#################
    #println("The given final points")
    #println(poly.final_config[1].x, ", ", poly.final_config[1].y, ", ", poly.final_config[1].z, ", ", poly.final_config[1].p)
    #Print the configurations and compare to what is calculated within error
    #println("The final points of the polynomial")
    #println(evaluate_poly(poly.x_coeffs, 0, poly.t), ", ", evaluate_poly(poly.y_coeffs, 0, poly.t), ", ", evaluate_poly(poly.z_coeffs, 0, poly.t), ", ", evaluate_poly(poly.p_coeffs, 0, poly.t))
    #initial velocities
    #println("The given final velocities")
    #println(poly.final_config[2].x, ", ", poly.final_config[2].y, ", ", poly.final_config[2].z, ", ", poly.final_config[2].p)
    #Print the configurations and compare to what is calculated within error
    #println("The final velocities of the polynomial")
    #println(evaluate_poly(poly.x_coeffs, 1, poly.t), ", ", evaluate_poly(poly.y_coeffs, 1, poly.t), ", ", evaluate_poly(poly.z_coeffs, 1, poly.t), ", ", evaluate_poly(poly.p_coeffs, 1, poly.t))
    #initial accelerations
    #println("The given final accelerations")
    #println(poly.final_config[3].x, ", ", poly.final_config[3].y, ", ", poly.final_config[3].z, ", ", poly.final_config[3].p)
    #Print the configurations and compare to what is calculated within error
    #println("The final accelerations of the polynomial")
    #println(evaluate_poly(poly.x_coeffs, 2, poly.t), ", ", evaluate_poly(poly.y_coeffs, 2, poly.t), ", ", evaluate_poly(poly.z_coeffs, 2, poly.t), ", ", evaluate_poly(poly.p_coeffs, 2, poly.t))
    ######################################################################333
    
    
    
    
    
    
    
    
    # Form a grid which shows the cells we check:
    N = round(Int64, ceil(get_grid_extent()/get_grid_resolution())); # Number of points in a dimension
    xy_check_grid = zeros(N,N); # 0 means not checked.
    xz_check_grid = zeros(N,N); # 0 means not checked.
    yz_check_grid = zeros(N,N); # 0 means not checked.
    
    #bool for checking if fine
    fine = true;
    for cell in poly.cells
        # Convert to subindices:
        subs = ind2sub((N,N,N), cell);
        #println(subs)
        if(subs[1] .> N || subs[2] .> N || subs[3].> N || subs[1]  .< 0 || subs[2]  .< 0 || subs[3]  .< 0)
            fine = false;
            break
        end
        #println("This is fine: ", fine)
        #println(subs)
        # Indicate that we check:
        xy_check_grid[subs[1],subs[2]] += 1;
        xz_check_grid[subs[1],subs[3]] += 1;
        yz_check_grid[subs[2],subs[3]] += 1;
    end

    #Only plot if possible

    if(fine)
        # plot the grid (floor level slice)
        figure(1,figsize=(6,6)); clf();    
        subplot(2,2,1);
        imshow(xy_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
        plot(xs+delta,width-ys-delta,color=:red);    
        title("XY plane");
        subplot(2,2,2);
        imshow(yz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
        plot(ys+delta,width-zs-delta,color=:red)
        title("YZ plane");
        subplot(2,2,3);
        imshow(xz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
        plot(xs+delta,width-zs-delta,color=:red)
        title("XZ plane");
        # plot again without the poly
        figure(2,figsize=(6,6)); clf();    
        subplot(2,2,1);
        imshow(xy_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
        #plot(xs+delta,width-ys-delta,color=:red);    
        title("XY plane");
        subplot(2,2,2);
        imshow(yz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
        #plot(ys+delta,width-zs-delta,color=:red)
        title("YZ plane");
        subplot(2,2,3);
        imshow(xz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
        #plot(xs+delta,width-zs-delta,color=:red)
        title("XZ plane");
        # Plot the polynomial:
        #plot(xs,ys,color=:red, linewidth=2);

    end
    
end


function test_tester()
    vel_lim = 1*2;
    accel_lim = 0.65/4*2;
    pt_lim = (get_grid_extent()*0.8);
    polyseg = connect_points([Point(rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1);Point((rand()-0.5)*vel_lim,(rand()-0.5)*vel_lim,(rand()-0.5)*vel_lim,(rand()-0.5)*vel_lim);Point((rand()-1)*accel_lim,(rand()-1)*accel_lim,(rand()-1)*accel_lim,(rand()-1)*accel_lim)], 
        [Point(rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1);Point((rand()-0.5)*vel_lim,(rand()-0.5)*vel_lim,(rand()-0.5)*vel_lim,(rand()-0.5)*vel_lim);Point((rand()-1)*accel_lim,(rand()-1)*accel_lim,(rand()-1)*accel_lim,(rand()-1)*accel_lim)], q_coeff, 5000.0, 100);

   polyseg = connect_points([Point(1.0,1.0,0.0,0.0);Point(1.0,1.0,0.0,0.0);Point(0.0,0.0,0.0,0.0)], 
        [Point(1.0,2.0,0.0,0.0);Point(0.0,0.0,0.0,0.0);Point(0.0,0.0,0.0,0.0)], q_coeff, 5000.0, 100);
    
    #polyseg = connect_points([Point(rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1);Point(rand()*vel_lim,rand()*vel_lim,rand()*vel_lim,rand()*vel_lim);Point(rand()*accel_lim,rand()*accel_lim,rand()*accel_lim,rand()*accel_lim)], 
        #[Point(rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1,rand()*pt_lim+1);Point(rand()*vel_lim,rand()*vel_lim,rand()*vel_lim,rand()*vel_lim);Point(rand()*accel_lim,rand()*accel_lim,rand()*accel_lim,rand()*accel_lim)], q_coeff);
    #polyseg.cells = collect(1:50:2000);
    println("Finished the polynomial connnection");
    println("Cost: ", polyseg.q)
    #println(polyseg.init_config)
    #println("Initial velocities")
    #println(polyseg.init_config[2])
    #println("Final velocities")
    #println(polyseg.final_config[2])
    check_poly_segment(polyseg)
end
