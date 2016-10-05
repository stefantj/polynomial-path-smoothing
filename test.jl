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

# Now we test:
function check_poly_segment(poly::poly_segment)
    # Form a grid which shows the cells we check:
    N = round(Int64, ceil(get_grid_extent()/get_grid_resolution())); # Number of points in a dimension
    xy_check_grid = zeros(N,N); # 0 means not checked.
    xz_check_grid = zeros(N,N); # 0 means not checked.
    yz_check_grid = zeros(N,N); # 0 means not checked.

    println("Cells: ", poly.cells);

    for cell in poly.cells
        # Convert to subindices:
        subs = ind2sub((N,N,N), cell);
        # Indicate that we check:
        xy_check_grid[subs[1],subs[2]] += 1;
        xz_check_grid[subs[1],subs[3]] += 1;
        yz_check_grid[subs[2],subs[3]] += 1;
    end

    # plot the grid (floor level slice)
    width = get_grid_extent();
    figure(1,figsize=(6,6)); clf(); 
    subplot(2,2,1);
    imshow(xy_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
    title("XY plane");
    subplot(2,2,2);
    imshow(yz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
    title("YZ plane");
    subplot(2,2,3);
    imshow(xz_check_grid[:,:]', cmap="gray", interpolation="none",extent=[0,width,0,width]);
    title("XZ plane");

    # Plot the polynomial:
    num_tsteps=100;
    times = linspace(0,poly.t,num_tsteps);
    xs = zeros(num_tsteps)    
    ys = zeros(num_tsteps);
    zs = zeros(num_tsteps);

    order = size(poly.x_coeffs,1);
println("Order is $order");
    for ord = 1:order
        xs += poly.x_coeffs[ord].*(times.^(ord-1))
        ys += poly.y_coeffs[ord].*(times.^(ord-1))
        zs += poly.z_coeffs[ord].*(times.^(ord-1))
    end
    #plot(xs,ys,color=:red, linewidth=2);

    # Try 3D plot
    figure(2); clf();
    plot3D(xs,ys,zs);
    xlim([0,width]);
    ylim([0,width]);
    zlim([0,width]);


    # Should add more checks such as confirming velocity
    #Print the velocities
end


function test_tester()
    polyseg = connect_points([Point(0.0,6.0,9.0,1.0);Point(0.0,-0.25,0.0,0.0);Point(0.0,0.0,0.0,0.0)], 
        [Point(0.0,2.0,5.0,0.0);Point(0.0,0.0,0.25,0.0);Point(0.0,0.0,0.0,0.0)], q_coeff);
    #polyseg.cells = collect(1:50:2000);
    println("out of connnect");
    check_poly_segment(polyseg)
end
