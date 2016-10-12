using PyPlot


# Used to store locations in the configuration space
type Point
   x::Float64
   y::Float64
   z::Float64
   p::Float64
end

type PolySol                  # For 1 segment per pair of points
    num_segs::Int64             # Number of segments
    times::Vector{Float64}      # Times for each segment
    x_coeffs::Array{Float64,1}  # X coefficients
    y_coeffs::Array{Float64,1}  # Y coefficients
    z_coeffs::Array{Float64,1}  # X coefficients
    p_coeffs::Array{Float64,1}  # Y coefficients
end




#Function evaluate_poly evaluates a given derivative of a polynomial at a certain value
#Assumptions
# One dimensional polynomial
#Inputs
# coeffs - a vector of the polynomial coefficients that form the polynomial from c1 + c2*t + c3*t^2 + ...
# deriv - the number representing the ith derivative of the polynomial, a value of 0 means the 0th derivative
# t - the value(s) at which the derivative is evaluated at
#Outputs
# x - the value of the polynomial
function evaluate_poly(coeffs::Vector{Float64}, deriv::Int64, t)
    #see how many coefficients there are
    p = length(coeffs);
    #initialize x to zero
    x = 0;
    #loop as many times are needed to go through the polynomial after it has been derived
    for e = 0:(p-deriv-1)
        #the looping variable will start at zero since the polynomial's first power will start with at 0 no matter what
        #create the coefficient indexer by adding one to the looping variable and andding the deriv number
        ci = e + deriv +1;
        #calculate the derivative caused constants by using factorials
        d = factorial(e + deriv) / factorial(e);
        #combine the terms through multiplication and add to x
        x += d*coeffs[ci]*t.^e
        #repeat the loop
    end
    #return the end value
    return x;
end




#Function verifyActuateablePath checks if the smooth path is feasible given robots limits
#Assumptions
# Robot follows the quadrotor model set out by Minimum  Snap  Trajectory  Generation  and  Control  for  Quadrotors - Mellinger
# Euler angles of Z-X-Y
# Polynomials of the same degree
#Inputs
# solution - an object containing points, and times
# max_vel - the maximum velocity that the robot is limited to in ros
# max_accel - the maximum position acceleration that the robot is limited to 
# max_jerk - the maximum position jerk that the robot is limited to
# max_motor_rpm - the maximum rpm that a motor can get
#Outputs
# did_pass - a boolean that is true when the path passes
# timeProbv - a vector of times where motion is infeasible based on velocity
# timeProbm - a vector of times where motion is infeasible based on velocity
#Function verifyActuateablePath checks if the smooth path is feasible given robots limits
#Assumptions
# Robot follows the quadrotor model set out by Minimum  Snap  Trajectory  Generation  and  Control  for  Quadrotors - Mellinger
# Euler angles of Z-X-Y
# Polynomials of the same degree
#Inputs
# solution - an object containing points, and times
# max_vel - the maximum velocity that the robot is limited to in ros
# max_accel - the maximum position acceleration that the robot is limited to 
# max_jerk - the maximum position jerk that the robot is limited to
# max_motor_rpm - the maximum rpm that a motor can get
#Outputs
# did_pass - a boolean that is true when the path passes
# timeProbv - a vector of times where motion is infeasible based on velocity
# timeProbm - a vector of times where motion is infeasible based on velocity
function verifyActuateablePath(solution::PolySol, max_vel::Float64, max_accel::Float64, max_jerk::Float64, max_motor_rpm::Float64)
    #Extract important information from the solution object
    degree = 2 + 2*(3-1) #create the degree of each polynomial assuming 2 pts for each
    num_poly = solution.num_segs;
    xcoeffs = solution.x_coeffs;
    ycoeffs = solution.y_coeffs;
    zcoeffs = solution.z_coeffs;
    pcoeffs = solution.p_coeffs;
    time_vec = solution.times;
    timeProbv = zeros(0,1);
    timeProbm = zeros(0,1);
    #Set some important constants for this function
    did_pass = true;
    time_res = 100 #The resolution segments 
    red_degree_by = 2; #The two is because the we are taking the derivative twice
    red_degree = degree - red_degree_by; #create a reduced degree to used for the dd and ddd calcs
    #Needed Constants for calculating the motor rpm to be on path
    z_w = [0,0,1]; #the up vector for the world frame
    Jxx = 0.0036404; #Values for the inertia matrix maybe can take from somewhere
    Jyy = 0.0053670;
    Jzz = 0.0030522;
    Jxy = 2.9264e-6;
    Jxz = 2.3411e-5;
    Jyz = 0.0001756;
    I = [Jxx Jxy Jxz
        Jxy Jyy Jyz
        Jxz Jyz Jzz];
    gravity = 9.8; 
    u2rpm_1 = 5096839.959225280;  #Values for an inverted allocation matrix to convert input into rpm
    u2rpm_2 = 51485858.53986801;  #Can probably get from code already coded
    u2rpm_3 = -51485858.53986803;
    u2rpm_4 = 330817430.2740964;
    u2rpms = [u2rpm_1 u2rpm_2 u2rpm_3 u2rpm_4 #inverted allocation matrix
        u2rpm_1 -u2rpm_2 -u2rpm_3 u2rpm_4
        u2rpm_1 u2rpm_2 -u2rpm_3 -u2rpm_4
        u2rpm_1 -u2rpm_2 u2rpm_3 -u2rpm_4];
    mass = 0.800; # in kilograms

    #####################
    #create variables to hold values initialize containers with zeros
    #Needed for z_B and u1
    xdd = zeros(num_poly*time_res);
    ydd = xdd;
    zdd = xdd;
    #Needed for speed check
    xd = xdd;
    yd = xdd;
    zd = xdd;
    #Needed for good plotting
    x = zeros(Float64,0,1);
    y = zeros(Float64,0,1);
    z = zeros(Float64,0,1);
    #Needed for a_dot
    xddd = xdd;
    yddd = xdd;
    zddd = xdd;
    #Needed for the inputs to the copter
    u1_vec = xdd;
    u2_vec = xdd;
    u3_vec = xdd;
    u4_vec = xdd;
    #Create a time vector for reporting purpose
    timeRep = zeros(Float64,0,1);
    #Needed initializations
    z_B = 0; #variable to hold the up body vector of copter
    #Needed for coefficients and angular acceleration
    yawdd = xdd;
    yaw = xdd;
    #create containers for problem points
    xprob = zeros(0,1);
    yprob = zeros(0,1);
    xddprob = zeros(0,1);
    yddprob = zeros(0,1);
    xdddprob = zeros(0,1);
    ydddprob = zeros(0,1);
    xddddprob = zeros(0,1);
    yddddprob = zeros(0,1);
    ######################Simple Method#####################################
    #Create a figure for debugging plots
    figure();
    #For each segment
    for seg = 1:num_poly
        #Create indexing range for the coefficients
        index_range = (1:degree)+degree*(seg-1);
        #Get the coefficients of the polynomial
        xcoeffs_s = xcoeffs[index_range];
        ycoeffs_s = ycoeffs[index_range];
        zcoeffs_s = zcoeffs[index_range];
        pcoeffs_s = pcoeffs[index_range];
        #Create the time vector for the polynomials, because of the way the problem is set up
        # every polynomial must be evaluated from 0 to the end of its range and then shifted
        t = collect(linspace(0,time_vec[seg+1]-time_vec[seg],time_res));
        #Generate a time vector for plotting and reporting
        timeRep = [timeRep; t+time_vec[seg]]
        #Calculate positions
        x = [x; evaluate_poly(xcoeffs_s,0,t)];
        y = [y; evaluate_poly(ycoeffs_s,0,t)];
        z = [z; evaluate_poly(zcoeffs_s,0,t)];
        #Calculate the velocities
        xd = evaluate_poly(xcoeffs_s,1,t);
        yd = evaluate_poly(ycoeffs_s,1,t);
        zd = evaluate_poly(zcoeffs_s,1,t);
        #Calculate the velocities that the robot will experience
        total_vel = sqrt(xd.^2 + yd.^2 + zd.^2);
        #Calculate the accelerations
        xdd = evaluate_poly(xcoeffs_s,2,t);
        ydd = evaluate_poly(ycoeffs_s,2,t);
        zdd = evaluate_poly(zcoeffs_s,2,t);
        yawdd = evaluate_poly(xcoeffs_s,2,t);
        #Calculate the total position acceleration
        total_accel = sqrt(xdd.^2 + ydd.^2 + zdd.^2);
        #Calculate the jerks
        xddd = evaluate_poly(xcoeffs_s,3,t);
        yddd = evaluate_poly(ycoeffs_s,3,t);
        zddd = evaluate_poly(zcoeffs_s,3,t);
        #Calculate the total position jerk
        total_jerk = sqrt(xddd.^2 + yddd.^2 + zddd.^2);
        #Calculate the snaps
        #Compare against the limits and store time and points
        #Velocity Limits
        timeProbv = [timeProbv; timeRep[find(total_vel .> max_vel)]];
        xprob = [xprob; evaluate_poly(xcoeffs_s,0,t[find(total_vel .> max_vel)])];
        yprob = [yprob; evaluate_poly(ycoeffs_s,0,t[find(total_vel .> max_vel)])];
        #Acceleration Limits
        timeProbv = [timeProbv; timeRep[find(total_accel .> max_accel)]];
        xddprob = [xddprob; evaluate_poly(xcoeffs_s,0,t[find(total_accel .> max_accel)])];
        yddprob = [yddprob; evaluate_poly(ycoeffs_s,0,t[find(total_accel .> max_accel)])];
        #Jerk Limits
        timeProbv = [timeProbv; timeRep[find(total_jerk .> max_jerk)]];
        xdddprob = [xdddprob; evaluate_poly(xcoeffs_s,0,t[find(total_jerk .> max_jerk)])];
        ydddprob = [ydddprob; evaluate_poly(ycoeffs_s,0,t[find(total_jerk .> max_jerk)])];
        #Note what time and position it occurs
        #Plot the graphs for debugging
        #figure();
        #Plot Problem points on position graph
        #plot(t+time_vec[seg],zddd);
        #title("Problem")
    end
    #Print debugging information
    #println(xprob)
    #Plot points where path is infeasible
    #if(!isempty(timeProbv))
    #    figure(); 
    #    xlabel("X"); ylabel("Y"); title("Portion of the Path that has Feasibility Problems");
    #    plot(x,y,color=:gray,linestyle=":"); #Path
    #    scatter(xprob, yprob); # Velocity Prob
    #    scatter(xddprob, yddprob, color=:green);
    #    scatter(xdddprob, ydddprob, color=:red);
    #    legend(["Path","Velocity", "Acceleration", "Jerk"], loc=0)
    #end
    #Return time locations where limits were exceeded.
    return timeProbv;
    ########################Complex method##################################
    #In the complex method the copter is only limited by the motor, arbitrary
    # limits can be applied should theoretically be unnecessary
    #For each segment
    #Calculate second derivatives of segments
    #Calculate z_B - body frame up vector
    #Calculate x_c
    #Calculate y_B
    #Calculate x_B
    #Calculate third derivatives
    #Calculate u1
    #Calculate h_w
    #Calculate yaws
    #Calculate w_bw
    #Calculate u2 u3 u4
    #Calculate the rpms
    #Compare to copter capabilities
    #Note what times and positions did the path become infeasible
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
#Local function to create Q
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

#Funtion occupancyCellChecker gets a path and finds the unique cell IDs in the occupancy grid that the path goes through.
#Assumptions
# An occupancy function exists that is called as follows ID = occupancy_get_id(x,y,z)
# Path starts at zero
# Solution has coefficient vectors of the same length
# The grid is [0, get_grid_extent] x [0, get_grid_extent]
# The dimension must be between 1 and 3 inclusive
# Cutting corners okay sometimes
#Inputs/Needed Variables
# xcf,ycf,zcf - coefficients of poly solution
# times - the times of the poly
# grid_resx - the resolution of the grid in the x direction
# grid_resy - the resolution of the grid in the y direction
# grid_resz - the resolution of the grid in the z direction
# dim - the dimension of the path to be checked
# timeStep - how much a time should be incremented to walk throught he plynomial path
# aggressParam - a parameter of how aggressive to check a path (0,infty), Closer to 0 means check every point on the path, Closer to infinity => check no points
#Outputs
# occupancy_vec - a vector of the occupancy IDs that the polynomial is characterized by
# outOfBounds - a boolean that is true if the polynomial goes outside of the grid
function occupancyCellChecker(xcf,ycf,zcf,times, grid_resx::Float64, grid_resy::Float64, grid_resz::Float64, dim::Int64, aggressParam::Float64, timeStep::Float64)
    #Check if the dimension is reasonable otherwise print error and exit;
    if(dim < 1 || dim > 3)
        println("Invalid dimension entered")
        return -1;
    end
    
    #Start outOfBounds as false
    outOfBounds = false;

    #Read in the coefficients into a matrix
    coeffMat = [xcf'; ycf'; zcf'];
    #Create a holder for derivatives, delta_t's, and pts when we get there in the for loop
    derivMat = zeros(dim);
    delta_t = zeros(dim);
    pts = zeros(3); #Three is used here since if dim is something other than 3 at this point the values are set to zero
    #prevPt = zeros(dim);
    occupancy_vec = zeros(Int64, 0,1);
    

    
    #Initialize/read in the resolution variables if needed
    #grid_resx = 0.05;
    #grid_resy = 0.05;
    #grid_resz = 0.05;
    #Initialize times assuming polynomials are solve at 0 time initial
    t = 0;
    timeFin = maximum(times) - minimum(times);
    #Time step will give you the resolution of the time steps in the loops
    #timeStep = timeFin/1000.0; ##############################################May want to make a variable
    #Calculate the distance time to travel in each direction before getting occupancy grid id using 1st order taylor series approximation
    dist_x = grid_resx * aggressParam;
    dist_y = grid_resy * aggressParam;
    dist_z = grid_resz * aggressParam;
    #Put into a matrix
    dist_to_travel = [dist_x; dist_y; dist_z]
    #Loop for 2 or 3 dimensions if dim >3
    for looper = 1:dim
        #Create previous point vector
        pts[looper] = evaluate_poly(coeffMat[looper,:],0,t);
        #Check if in bounds
        if(pts[looper] < 0 || pts[looper] > get_grid_extent())
            outOfBounds = true;
        end
        #Create a variable to change temporarily
        t_new = t;
        counter = 0;
        #Calculate the time at which the change in distance is more than our distance to travel
        while(abs(pts[looper] - evaluate_poly(coeffMat[looper,:],0,t_new)) < dist_to_travel[looper] && counter < 1000)
            #If we haven't gotten higher than what we want to travel increment the time
            t_new += timeStep;
            #println(t_new)
            counter += 1;
        end
        #println("past first while")
        #Set the new delta_t
        delta_t[looper] = abs(t_new - t);
        #evaluate the polynomial first direvative at current time
        #derivMat[looper] = evaluate_poly(coeffMat[looper, :], 1, t);
        #if(derivMat[looper] == 0)
        #    derivMat[looper] = 1;
        #end
        #looper = 1 -> x, 2 -> y, 3 -> z
        #find the required delta_t's needed to go the set distance before checking occupancy
        #delta_t[looper] = abs(dist_to_travel[looper] / derivMat[looper]);
    end

    #Create a counter
    counter1 =1;

    #Check if at the or past the end time and loop if not passed it
    while(t < timeFin ) #&& counter1 < 100)
        #Find the smallest one
        changeDelta, deltaIndex = findmin(delta_t)
        #Move the time by the smallest one
        t += changeDelta;
        #Decrease the delta t's for the others
        for looper = 1:dim
            #Subract the smallest delta t from the other deltas
            delta_t[looper] -= changeDelta;
        end
        #Find the x,y,z of the current time
        for looper = 1:dim
            pts[looper] = evaluate_poly(coeffMat[looper, :], 0, t);
            #Check if in bounds
            if(pts[looper] < 0 || pts[looper] > get_grid_extent())
                outOfBounds = true;
            end
        end
        #Find occupancy ID at current point by adding to a vector
        occupancy_vec = [occupancy_vec; occupancy_get_id(pts[1], pts[2], pts[3])];
        #println(occupancy_vec)
        #println(changeDelta)
        #Evaluate the new delta_t for the dimension
        for looper = 1:dim
            if(delta_t[looper] == 0)
                #Create previous point vector
                #pts[looper] = evaluate_poly(coeffMat[looper,:],0,t);
                #Create a variable to change temporarily
                t_new = t;
                counter = 0;
                #Calculate the time at which the change in distance is more than our distance to travel
                while(abs(pts[looper] - evaluate_poly(coeffMat[looper,:],0,t_new)) < dist_to_travel[looper] && counter < 1000)
                    #If we haven't gotten higher than what we want to travel increment the time
                    t_new += timeStep;
                    counter += 1;
                end
                #Set the new delta_t
                delta_t[looper] = abs(t_new - t);
            end
        end
        #derivMat[deltaIndex] = evaluate_poly(coeffMat[deltaIndex, :], 1, t);
        #Make sure the derivative is not zero
        #if(derivMat[deltaIndex] == 0)
        #    derivMat[deltaIndex] = 1;
        #end
        #delta_t[deltaIndex] = abs(dist_to_travel[deltaIndex]/derivMat[deltaIndex]);
        #Increment counter
        #counter1 += 1;
    end

        
    #Check the final point just in case based on passed dimension
    #println("final point")
    if(dim == 3)
        occupancy_vec = [occupancy_vec; occupancy_get_id(evaluate_poly(coeffMat[1,:],0,timeFin),
            evaluate_poly(coeffMat[2,:],0,timeFin),evaluate_poly(coeffMat[3,:],0,timeFin))];
    elseif(dim == 2)
        occupancy_vec = [occupancy_vec; occupancy_get_id(evaluate_poly(coeffMat[1,:],0,timeFin),
            evaluate_poly(coeffMat[2,:],0,timeFin),0)];
    else
        occupancy_vec = [occupancy_vec; occupancy_get_id(evaluate_poly(coeffMat[1,:],0,timeFin),0,0)];
    end
    
    #Return the vector of unique occupancy grids
    return unique(occupancy_vec), outOfBounds;
end

# Function which returns the cell id of the point (xyz). Returns an integer.
function occupancy_get_id(x,y,z)
    width = get_grid_extent();
    res   = get_grid_resolution();
    n = round(Int64,ceil(width/res));
    #println(x,y,z)
    return sub2ind((n,n,n),floor(Int64,x/res)+1,floor(Int64,y/res)+1,floor(Int64,z/res)+1)
end

function get_grid_extent()
    return 10;
end
function get_grid_resolution()
    return 0.1
end