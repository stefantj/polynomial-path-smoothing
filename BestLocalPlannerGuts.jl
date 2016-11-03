## constructConstr
#Description - Create constraint matrices that follow convention and the desired order.
#Assumptions
# Every Polynomial will be made the same
# All soft and free constraints will be end derivatives
# Constraints ordered in fixed, free, soft order with pos, vel,.. and init then final
#Inputs
# prob - a PathProblem type that contains the following important values:
#   start_config::Array{Float64,2}   # Initial conditions with pos, vel, acc, jerk, ...
#   end_config::Array{Float64,2}     # Final conditions with pos, vel, acc, jerk,... in that order
#   soft_constr::Array{Bool,1}       # A vector that determines whether the end configurations are soft or not
#   PconstrFixed::Array{Float64,2}   # All fixed constraints
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   PtimeIndex::Array{Int64,1}       # The vector of the times at which the constraints apply
#   PconstraintOrders::Array{Int64,1}# The order corresponding to each constraint in the total constriant vector
#   PconstraintSoft::Array{Float64,2}# Soft constraints ordered by ascending order
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
#   dof::Int64                       # The extra degrees that would be completely free to optimize over
#   Pdegree::Int64                   # The degree of the polynomial to be made
#Outputs
# Updates prob's following members
#   PconstrFixed::Array{Float64,2}   # All fixed constraints
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   PtimeIndex::Array{Int64,1}       # The vector of the times at which the constraints apply
#   PconstraintOrders::Array{Int64,1}# The order corresponding to each constraint in the total constriant vector
#   PconstraintSoft::Array{Float64,2}# Soft constraints ordered by ascending order
#   Pdegree::Int64                   # The degree of the polynomial to be made
function constructConstr(prob::PathProblem)
    #Determine the total degrees of freedom needed for the poly
    numFinal = size(prob.end_config,2);      #number of final conditions
    numInit = size(prob.start_config,2);     #number of initial conditions
    degree = numInit + numFinal + prob.dof;
    #Create a vector of indices for where softs do and do not occur
    softIndex = find(prob.soft_constr .== true);
    notSoftIndex = find(prob.soft_constr .== false);
    #Determine the soft constraints and save them
    constrSoft = zeros(2+prob.isDim3,length(softIndex));
    for polys =1:(2+prob.isDim3)
        #Create constrSoft in order
        constrSoft[polys,:] = prob.end_config[polys,softIndex];
    end
    #Subtract dofExtra and the number of soft constraints from the total degree to get the size of constrFixed
    #Create the fixed constraint vector

## STEFAN: Changed some dimensions to what I think looks right. Double check though.
    constrFixed =  zeros(2+prob.isDim3,degree-length(softIndex)-prob.dof);
    for i = 1:(2+prob.isDim3)
        if(size(notSoftIndex,1) > 0)
            constrFixed[i,:] = [prob.start_config[i,:] prob.end_config[i,notSoftIndex]];
        else
            constrFixed[i,:] = [prob.start_config[i,:]];
        end
    end
    #Create the free vector with soft constraints at the end
    constrFree =  zeros(2+prob.isDim3,length(softIndex)+prob.dof);
    for i = 1:(2+prob.isDim3)
        constrFree[i,:] = [zeros(1,prob.dof) constrSoft[i,:]]
    end
    #Create the time index vector with assumption that all free variables are end derivatives
    timeIndex = ones(degree);
    timeIndex[1:numInit] = zeros(numInit);  
    #Create the order vector assuming that the index of the original constraint is one more than its actual order
    constrOrder = [
        collect((1:numInit)-1)                        #init fixed
        notSoftIndex-1                                #final fixed
        (1:prob.dof)+size(prob.soft_constr,1)-1       #additional degrees of freedom
        softIndex-1                                   #final soft
    ]
    #Update the path problem
    prob.PconstrFixed = constrFixed;
    prob.PconstrFree = constrFree;
    prob.PtimeIndex = timeIndex;
    prob.PconstraintOrders = constrOrder;
    prob.Pdegree = degree;
    prob.PconstraintSoft = constrSoft;
    return prob;
end


## addDirectedSpeed!
#Description - Adds a directed velocity vector to the end configuration of the robot
#Assumptions
# Want to go as fast as possible 
# Not close to goal point
#Inputs
# prob - a PathProblem type that contains the following important values:
#   end_config::Array{Float64,2}     # Final conditions with pos, vel, acc, jerk,... in that order
#   soft_constr::Array{Bool,1}       # A vector that determines whether the end configurations are soft or not
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
#Outputs
# Updates prob's following members
#   end_config::Array{Float64,2}     # Final conditions with pos, vel, acc, jerk,... in that order
#   soft_constr::Array{Bool,1}       # A vector that determines whether the end configurations are soft or not
function addDirectedSpeed!(prob::PathProblem, max_vel)
    #Initialize a max_vel_vec
    max_vel_vec = [];
    #In whichever D (2 or 3)
    for i = 1:(2+prob.isDim3)
        max_vel_vec = [max_vel_vec; (prob.end_config[i,1] - prob.start_config[i,1])];
    end
    #Make sure that if max_vel_vec is zero you do not get NAN
    if(max_vel_vec != zeros(2+prob.isDim3))
        max_vel_vec = normalize!(max_vel_vec)*max_vel
    end
    #Add as a soft constraint being careful about what exists
    if(size(prob.end_config,2) <= 1)
        prob.end_config = [prob.end_config[1:(2+prob.isDim3),1] max_vel_vec];
        prob.soft_constr = [prob.soft_constr[1]; true];
    else
        prob.end_config[1:(2+prob.isDim3),2] = max_vel_vec;
        #prob.soft_constr[2] = true;
    end
end

## constr_Ainvs
#Description - a function for deriving the entire A corresponding to constraint orders at given times for poly of given degree
# on a row by row basis. A is then inverted and output
#Assumptions
# Order and time are of the length degree
#Inputs
# order - a vector of the orders of constraints to construct A
# time - a vector of the times at which the constraints take place
# degree - the degree of the polynomial to be solved
#Outputs
# inv(A) - the inverse of the A matrix formed by the inputs
function constr_Ainv(order, time, degree)
    #Initialize an A matrix
    A = zeros(degree,degree);
    for i = 1:degree
        #Construct a row by row with legacy code that works well
        A[i,:] = constr_order(order[i], time[i], degree);
    end
    return inv(A);
end

## constr_order
#Description - Function for deriving the row of A corresponding to constraint of order at time for poly of degree.
#Assumptions
# Order and time are of the length degree
#Inputs
# order - an order number of constraint
# time - the time at which the constraint take place
# degree - the degree of the polynomial to be solved
#Outputs
# row_vec - the row_vec of A corresponding to the inputs
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

## form_Q
#Description - Function for making Jacobian for free derivative variable optimization. See Richter for more information.
#Assumptions
# Q_coeffs is the same size as the degree of the polynomial to be optimized
#Inputs
# Q_coeffs - the weights on the derivative terms
# t - the total time of the polynomial
#Outputs
# Q_mat - the Jacobian for the derivative 
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

## checkQcoeffs
#Description - makes sure that the q_coeff vector is the same size as the problem's degree. Additions are only zeros.
#Assumptions
# q_coeff is an array of ints
#Inputs
# q_coeff - a vector of waits for the Jacobian of derivative constraint costs
# degree - the degree that the q_coeff must have the size of
#Outputs
# an updated q_coeff of same size as degree
function checkQcoeffs(q_coeff, degree)
    for i = 1:abs(size(q_coeff,1)-degree)
        if(size(q_coeff,1)-degree < 0)
            #Add a 0 if q_coeff is too small
            q_coeff = [q_coeff; 0.0]
        elseif(size(q_coeff,1)-degree > 0)
            #Subtract if too big
            q_coeff = q_coeff[1:end-1];
            #Else nothing should happen
        end
    end
    return q_coeff;
end

## form_OptimizeMat
#TODO - make this function a void that directly changes solvStuff
#Description - forms the optimizeMatrix that translates fixed constraints into the best free constraint choices for minimizing the
# derivactive costs.
#Assumptions
# Only one segment
#Inputs
# prob - the PathProblem object with the following values:
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   PconstraintSoft::Array{Float64,2}# Soft constraints ordered by ascending order
# tuning - the TuningParams object with the following values:
#   softConstrWeights::Array{Float64,1}# The vector of weights to make soft constraints harder
# solvStuff - a PolyPathSolver object with the following values:
#   PA_inv::Array{Float64,2}           # The matrix that transforms constraints to coefficients A_inv*d=p
#   PQ::Array{Float64,2}               # Matrix representation of the q_coeff weights
#Outputs
# optimizeMat - the matrix that transforms the fixed constraints to optimal values of free constraints w.r.t. derivative costs
function form_OptimizeMat(prob::PathProblem,tuning::TuningParams,solvStuff::PolyPathSolver)
    #Create the R matrix for the poly 
    R = solvStuff.PA_inv'*(solvStuff.PQ*solvStuff.PA_inv);
    #Create a holder for the number of free and soft constraints
    num_free = size(prob.PconstrFree,2);
    num_soft = size(prob.PconstraintSoft,2);
    #Create the soft constraint weight matrix
    softWeight = zeros(num_free,num_free)
    #TODO: make tuning a pointer so that this next change doesn't have to constantly be redone
    #Make sure that there are enough weights
    tuning.softConstrWeights = checkQcoeffs(tuning.softConstrWeights, num_soft)
    #Add the weights from the bottom right to upper left
    for i = (1:num_soft)-1
        softWeight[end-i,end-i] = tuning.softConstrWeights[end-i];
    end
    #Create the RppInv term in the analytical optimization
    RppInv = inv(R[(end-num_free+1):end, (end-num_free+1):end]+softWeight)
    optimizeMat = -RppInv*R[1:(end-num_free),(end-num_free+1):end]';
    return optimizeMat
end

## updateFreeConstr
#Description - will do the matrix multiplication of optimizeMat and the fixed constraints to get the optimal free constraints
# w.r.t. the derivative costs.
#Assumptions
# All the dimensions of the problem have the same number of soft constraints and free constraints
#Inputs
# prob - the PathProblem object with the following values:
#   PconstrFixed::Array{Float64,2}   # All fixed constraints
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   PconstraintSoft::Array{Float64,2}# Soft constraints ordered by ascending order
# solvStuff - a PolyPathSolver object with the following values:
#   PoptimizeMatrix::Array{Float64,2}  # Gives the optimal free/soft constraints given fixed constraints
#Outputs
# PconstrFree - the updated free constraints
function updateFreeConstr(prob::PathProblem,solvStuff::PolyPathSolver)
    #Create a temporary containter to hold the soft constraints to add
    addSofts = zeros(size(prob.PconstrFree));
    #For each dimension
    for i = 1:(2+prob.isDim3)
        #Non soft constraints will have zeros
        addSofts[i,end-size(prob.PconstraintSoft,2)+1:end] = prob.PconstraintSoft[i,:];
        #The actual calculation is Opt_mat * d_f + d_softs
### STEFAN: Changed dimensions to what I think looks right. Double check 
        prob.PconstrFree[i,:] = prob.PconstrFixed[i,:]*solvStuff.PoptimizeMatrix'  + addSofts[i,:];
    end
    return prob.PconstrFree;
end


## solvePolysInitially
#Description - will do the matrix multiplication of optimizeMat and the fixed constraints to get the optimal free constraints
# w.r.t. the derivative costs.
#Assumptions
# All the dimensions of the problem have the same number of soft constraints and free constraints
#Inputs
# prob1 - the PathProblem object with the following values:
#   PconstrFixed::Array{Float64,2}   # All fixed constraints
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3  
#   Pdegree::Int64                   # The degree of the polynomial to be made
# solvStuff - a PolyPathSolver object with the following values:
#   PA_inv::Array{Float64,2}         # The matrix that transforms constraints to coefficients A_inv*d=p
#Outputs
# coeffs - the polynomial coefficients in a matrix
function solvePolysInitially(prob1,solvStuff)
    #initialize the coeeffs variable
    coeffs = zeros(2+prob1.isDim3,prob1.Pdegree)
    #for each dimension
    for i = 1:(2+prob1.isDim3)
println(size(solvStuff.PA_inv), size(prob1.PconstrFixed[i,:]'), size(prob1.PconstrFree[i,:]'));
### STEFAN: Changed dimensions to what I think looks right. Double check 
        #p = A^-1 * d
        coeffs[i,:] = (solvStuff.PA_inv * [prob1.PconstrFixed[i,:]'; prob1.PconstrFree[i,:]'])';
    end
    return coeffs;
end

## createDerivMat
#Description - will create a matrix that will map a polynomial's coefficient values to the values of that polynomial's user specified
# derivative.
#Assumptions
# Coefficients are order p1t^0 + p2*t^1+.....
#Inputs
# deg - the nnumber of coefficients of the polynomial
# deriv - the derivative to make the matrix for 1 = 1st, 2 = 2nd
#Outputs
# mat - the matrix that ill map a polynomial's coefficient values to the values of that polynomial's derivative
function createDerivMat(deg::Int64, deriv::Int64)
    #initialize matrix to zero
    mat = zeros(deg,deg);
    #loop as many times are needed to go through the polynomial after it has been derived
    for e = 0:(deg-deriv-1)
        #calculate the derivative caused constants by using factorials
        derivConst = factorial(e + deriv) / factorial(e);
        #insert the derivative caused constants at the appropriate diagonal position
        mat[e+deriv+1, e+deriv+1] = derivConst;
        #TODO: make better
    end
    return mat;
end

## occupancyCellChecker 
#Description gets a path and finds the unique cell IDs in the occupancy grid that the path goes through.
#Assumptions
# An occupancy function exists that is called as follows ID = occupancy_get_id(x,y,z)
# Path starts at zero time
# Solution has coefficient vectors of the same length
# The grid is [0, get_grid_extent] x [0, get_grid_extent]
# There is uniform resolution in each direction
# The dimension must is either 2 or 3
# Cutting corners okay sometimes
#Inputs/Needed Variables
# sol - a PathSol object that has the following values
#   totTime::Float64            # Total time of the polynomial path in seconds
#   coeffs::Array{Float64,2}    # X, Y, Z, and Yaw coefficients
# prob - a PathProblem object with the following values
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
#   grid_extent::Float64             # The extent of the cost map in meters
#   grid_resolution::Float64         # The resolution of the cost map in meters
# tuning - a TuningParams object that has the following variables
#   timeStep::Float64                  # The increment in time to determine the cells that a path goes through
#   aggressParam::Float64              # Closer to 0 means check every point on the path, infinity => check no points
#Outputs
# occupancy_vec - a vector of the occupancy IDs that the polynomial is characterized by (not unique but can be)
# outOfBounds - a boolean that is true if the polynomial goes outside of the grid
function occupancyCellChecker(sol::PathSol, prob::PathProblem, tuning::TuningParams)
   
    #Create the dim variable
    dim = 2+prob.isDim3;
    #Start outOfBounds as false
    outOfBounds = false;

    #Read in the coefficients into a matrix
    coeffMat = sol.coeffs
    #Create a holder for delta_t's and pts when we get there in the for loop
    delta_t = zeros(dim);
    pts = zeros(3); #Three is used here since if dim is something other than 3 at this point the values are set to zero
    occupancy_vec = zeros(Int64, 0,1);
   
    #TODO: Make this a do while
    #Initialize times assuming polynomials are solve at 0 time initial
    t = 0;
    timeFin = sol.totTime;
    #Time step will give you the resolution of the time steps in the loops
    timeStep = tuning.timeStep;
    #Calculate the distance time to travel in a direction before getting occupancy grid id 
    dist = prob.grid_resolution * tuning.aggressParam;
    #Loop for 2 or 3 dimensions
    for looper = 1:dim
        #Create previous point vector
### STEFAN: Changed dimensions to what I think looks right. Double check 
        pts[looper] = evaluate_poly(vec(coeffMat[looper,:]),0,t);
        #Create an initial dist that should be checked the first time 
        distInit = minimum( [mod(pts[looper],prob.grid_resolution),
                             prob.grid_resolution-mod(pts[looper],prob.grid_resolution)])
        #Check if in bounds and stop the function early
        if(pts[looper] < 0 || pts[looper] > prob.grid_extent)
            outOfBounds = true;
            return unique(occupancy_vec), outOfBounds;
        end
        #Create a variable to change temporarily
        t_new = t;
        counter = 0;
        #Calculate the time at which the change in distance is more than our distance to travel
        while(abs(pts[looper] - evaluate_poly(vec(coeffMat[looper,:]),0,t_new)) < distInit && counter < 1/timeStep)
            #If we haven't gotten higher than what we want to travel increment the time
            t_new += timeStep;
            #println(t_new)
            counter += 1;
        end
        #Set the new delta_t
        delta_t[looper] = abs(t_new - t);
    end

    #Check if at the or past the end time and loop if not 
    while(t < timeFin) 
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
            pts[looper] = evaluate_poly(vec(coeffMat[looper, :]), 0, t);
            #Check if in bounds if not leave the function early
            if(pts[looper] < 0 || pts[looper] > prob.grid_extent)
                outOfBounds = true;
                return unique(occupancy_vec), outOfBounds;
            end
        end
        #Find occupancy ID at current point by adding to a vector
        occupancy_vec = [occupancy_vec; occupancy_get_id(pts[1], pts[2], pts[3], prob)];
        #Evaluate the new delta_t for the dimension
        for looper = 1:dim
            if(delta_t[looper] == 0)
                #Create a variable to change temporarily
                t_new = t;
                counter = 0;
                #Calculate the time at which the change in distance is more than our distance to travel
                while(abs(pts[looper] - evaluate_poly(vec(coeffMat[looper,:]),0,t_new)) < dist && counter < 1000)
                    #If we haven't gotten higher than what we want to travel increment the time
                    t_new += timeStep;
                    counter += 1;
                end
                #Set the new delta_t
                delta_t[looper] = abs(t_new - t);
            end
        end
    end
    #Check the final point just in case based on passed dimension
    #println("final point")

    if(dim == 3)
        occupancy_vec = [occupancy_vec; occupancy_get_id(vec(evaluate_poly(coeffMat[1,:]),0,timeFin),
                                                         evaluate_poly(vec(coeffMat[2,:]),0,timeFin),
                                                         evaluate_poly(vec(coeffMat[3,:]),0,timeFin), prob)];
    elseif(dim == 2)
        occupancy_vec = [occupancy_vec; occupancy_get_id(evaluate_poly(vec(coeffMat[1,:]),0,timeFin),
                                                         evaluate_poly(vec(coeffMat[2,:]),0,timeFin),
                                                         0,prob)];
    end
    #Return the vector of unique occupancy grids (could be non unique)
    return unique(occupancy_vec), outOfBounds;
end

## occupancy_get_id
#Description - returns the cell id of the point (xyz). Returns an integer.
#Assumptions
# All point are within the map
#Inputs
# x,y,z - the respective points on a path
# problem1 - a PathProblem object with the following members:
#   grid_extent::Float64             # The extent of the cost map in meters
#   grid_resolution::Float64         # The resolution of the cost map in meters
#Outputs
# The integer index to a 3d array
function occupancy_get_id(x,y,z, problem1)
    width = problem1.grid_extent;
    res   = problem1.grid_resolution;
    n = round(Int64,ceil(width/res));
    tempx = x
    tempy = y
    tempz = z
    return sub2ind((n,n,n),floor(Int64,tempx/res)+1,floor(Int64,tempy/res)+1,floor(Int64,tempz/res)+1)
end

## simpleVerifyFeas
#Description - checks if the smooth path is feasible given robots limits on speed, accel, and jerk
#Assumptions
# Polynomials of the same degree
# Only one segment
#Inputs
# sol - a PathSol object containing:
#   totTime::Float64            # Total time of the polynomial path in seconds
#   coeffs::Array{Float64,2}    # X, Y, Z, and Yaw coefficients
# tuning - a TuningParams object that contains: 
#   max_vel::Float64                   # Max total velocity the path is restricted to
#   max_accel::Float64                 # Max total acceleration the path is restricted to
#   max_jerk::Float64                  # Max total jerk the path is restricted to
#   timeRes::Int64                     # How many points for determining limits and constraints
#Outputs
# timeProbv - a vector of values that did not pass
# timeProbm - a vector of numbers 1-3 corresponding to where motion is infeasible based on 1 velocity, 2 accel, and 3 jerk
#Function verifyActuateablePath checks if the smooth path is feasible given robots limits
function simpleVerifyFeas(sol::PathSol, tuning::TuningParams)
    #Create empty holders at first
    timeProbv = zeros(0,1);
    timeProbm = zeros(0,1);
    #Create a matrix of the limits
    max_lims = [tuning.max_vel;
                tuning.max_accel;
                tuning.max_jerk]
    #create a dim variable
    dim = size(sol.coeffs,1);

    #Create the time vector 
    t = collect(linspace(0,sol.totTime,tuning.timeRes));
    
    #Initialize a container for checking all the limits
    checkingMat = zeros(dim,length(t))
    for j = 1:3
        #Calculate the derivative
        #initialize the summer
        totalSum = 0;
        for i = 1:dim
            #calculate the respective derivative
            checkingMat[i,:] = evaluate_poly(vec(sol.coeffs[i,:]),j,t);
            #sum its square
            totalSum += checkingMat[i,:].^2;
        end
        #Sqaure root the sum and compare
        totalSum = sqrt(totalSum);

        timeProbv = [timeProbv; totalSum[find(totalSum .> max_lims[j]+tuning.precisionVel)]];
        timeProbm = [timeProbm; 0*find(totalSum .> max_lims[j]+tuning.precisionVel)+j]
    end
    return timeProbv, timeProbm;
end

## initializeTime
#Description - sets the time to be as fast as possible given the distance to travel and the max_vel. This is done to avoid an
# observed phenomenon of the solution shooting out of bounds if the initial time lead to huge velocities, accels, and jerks.
#Assumptions
# This will avoid the phenomenon (not sure this is valid for all cases especially if robot is going in the wrong direction to start)
#Inputs
# prob - a PathProblem object that has:
#   start_config::Array{Float64,2}   # Initial conditions with pos, vel, acc, jerk, ...
#   end_config::Array{Float64,2}     # Final conditions with pos, vel, acc, jerk,... in that order
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
# tuning - a TuningParams object that has:
#   timeIncrease::Float64              # The increment by which to increase time if path is found infeasible
#   max_vel::Float64                   # Max total velocity the path is restricted to
#Outputs
# totTime - the reasonable total time to start the solution off at
function initializeTime(prob, tuning)
    #Initialize
    totTime = 0.0;
    for i = 1:(2+prob.isDim3)
        #Sum sqaure distances
        totTime += (prob.start_config[i,1]-prob.end_config[i,1])^2;
    end
    #sqrt and then divide by max_vel and then subtract the timeIncrease that will be added immediately after this
    totTime = sqrt(totTime)
    totTime /= tuningI.max_vel;
    totTime -= tuningI.timeIncrease;
    return totTime;
end

## debugPlot
#Description - plots the path and various derivative magnitudes with no background and can be used even if out of bounds
#Inputs
# sol - PathSol object with:
#   totTime::Float64            # Total time of the polynomial path in seconds
#   coeffs::Array{Float64,2}    # X, Y, Z, and Yaw coefficients
# prob - PathProblem object with:
#   isDim3::Bool                # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
# tuning - TuningParams object with:
#   timeRes::Int64              # How many points for determining limits and constraints
#Outputs
# Plots... with labels son! Aight! Cool...
function debugPlot(sol, prob, tuning)
    #Create times according to resolution
    plotTimes = linspace(0,sol.totTime,tuning.timeRes)
    #First figure is position
    figure(1)
    if(prob.isDim3)
        plot3D(   evaluate_poly(sol.coeffs[1,:],0,plotTimes), #x
                evaluate_poly(sol.coeffs[2,:],0,plotTimes), #y
                evaluate_poly(sol.coeffs[3,:],0,plotTimes)); #z
        #Extra label for 3D: rad!
        zlabel("Z (m)")
    else
        plot(   evaluate_poly(vec(sol.coeffs[1,:]),0,plotTimes), #x
                evaluate_poly(vec(sol.coeffs[2,:]),0,plotTimes)) #y 
    end
    #Labels Homeboi!
    title("Position")
    xlabel("X (m)")
    ylabel("Y (m)")
    #Titles and ylabels for reference in the for loop boi!
    titles = ["Velocity Magnitude", "Acceleration Magnitude", "Jerk Magnitude"]
    ylabels = ["Velocity Magnitude (m/s)", "Acceleration Magnitude(m/s^2)", "Jerk Magnitude (m/s^3)"]
    for i =1:3
        figure(i+1)
        if(prob.isDim3)
            #Various derivatives are sqaured summed and sqrt'ed
            plot(   plotTimes, sqrt(evaluate_poly(sol.coeffs[1,:],i,plotTimes).^2 + #x
                                    evaluate_poly(sol.coeffs[2,:],i,plotTimes).^2 + #y
                                    evaluate_poly(sol.coeffs[3,:],i,plotTimes).^2)) #z
        else
            plot(   plotTimes, sqrt(evaluate_poly(sol.coeffs[1,:],i,plotTimes).^2 + #x
                                    evaluate_poly(sol.coeffs[2,:],i,plotTimes).^2)) #y
        end
        title(titles[i])
        ylabel(ylabels[i])
        #Time label is same for all plots
        xlabel("Time (s)")
    end
end

## costFunc
#Description - the functionalized form of the objective function to be minimized with Optim.jl
#Assumptions
# Soft constraints are at the end
#Inputs
# dF - the free derivatives that can be varied in the optimization this must be separated for the optimization to work
#  the order of df is pressumed to be x derivs low to high, then y, and z if the dimension is three
# sol - the PathSol object with:
#   totTime::Float64                   # Total time of the polynomial path in seconds
#   coeffs::Array{Float64,2}           # X, Y, Z, and Yaw coefficients
#   cells::Array{Int64,1}              # The indices of the costmap that the path goes through
# prob - a PathProblem object that has:
#   PconstrFixed::Array{Float64,2}     # All fixed constraints
#   PconstrFree::Array{Float64,2}      # A holder for all soft and free contraints
#   isDim3::Bool                       # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
#   Pdegree::Int64                     # The degree of the polynomial to be made
# solvStuff - a PolyPathSolver object with: 
#   PA_inv::Array{Float64,2}           # The matrix that transforms constraints to coefficients A_inv*d=p
#   PQ::Array{Float64,2}               # Matrix representation of the q_coeff weights
# tuning - the TuningParams object with
#   softConstrWeights::Array{Float64,1}# The vector of weights to make soft constraints harder
#   numericVelWeight::Float64          # The weight for the cost for not being at max_vel at any time
#   obstacleWeight::Float64            # The weight for the cost of going through a cell in the costmap
#   derivativeWeight::Float64          # The weight for the combined derivative cost term
#   timeWeight::Float64                # Weight applied to time
#   timeRes::Int64                     # How many points for determining limits and constraints
#Outputs
# cost - the cost of the given path
function costFunc(dF, sol::PathSol, prob::PathProblem, solvStuff::PolyPathSolver, tuning::TuningParams)
    #Create a dim variable
    dim = 2+prob.isDim3;

    #Create a combined vector of fixed and free constraints and the resulting poly coefficients
    d = zeros(dim,prob.Pdegree)
    p = zeros(dim,prob.Pdegree)
### STEFAN: Changed dimensions to what I think looks right. Double check 
    for i = 1:dim
        d[i,:] = [prob.PconstrFixed[i,:] dF[(1:size(prob.PconstrFree,2))+(i-1)*size(prob.PconstrFree,2)]'];
        p[i,:] = solvStuff.PA_inv * (d[i,:]');
    end
    #Create the derivative cost
    derivCost = 0;
    for i = 1:dim
        derivCost += (  p[i,:]) *      # p^T
                        solvStuff.PQ *  # Q
                       (p[i,:]');        # p
    end
    #Divide by the number of derivative terms and then apply the weight
    #Note with a weight of 1 the cost is in the 100s
    derivCost = derivCost*tuning.derivativeWeight / dim #there may be an issue with dividing by and integer

    #Create the velocity cost
    #Create a time vector 
    timesCheck = collect(linspace(0,sol.totTime,tuning.timeRes));
    #Sum the square errors in velocity magnitudes
    if(prob.isDim3)
        velCost = sum((sqrt(evaluate_poly(vec(p[1,:]),1,timesCheck).^2 + # x
                            evaluate_poly(vec(p[2,:]),1,timesCheck ).^2 + # y
                            evaluate_poly(vec(p[3,:]),1,timesCheck ).^2)- # z
                       tuning.max_vel).^2) ;                       # minus the max_vel and square
    else
        velCost = sum((sqrt(evaluate_poly(vec(p[1,:]),1,timesCheck).^2 + # x
                            evaluate_poly(vec(p[2,:]),1,timesCheck).^2)- # y
                       tuning.max_vel).^2) ;                       # minus the max_vel and square
    end
    #Add a weight
    #Note that Vel cost is usually around order of 100 with a weight of 1
    velCost *= tuning.numericVelWeight;
    
    #Add an accel cost for going above the limit
    #Sum the square errors in acceleration magnitudes
    if(prob.isDim3)
        accelCostTemp = sqrt(evaluate_poly(p[1,:],2,timesCheck).^2 + evaluate_poly(p[2,:],2,timesCheck).^2 + evaluate_poly(p[3,:],2,timesCheck).^2)
    else
        accelCostTemp = sqrt(evaluate_poly(p[1,:],2,timesCheck).^2 + evaluate_poly(p[2,:],2,timesCheck).^2)
    end
    #Add a cost for being above and a no cost for being below
    accelCost = sum((accelCostTemp -  tuning.max_accel*tuning.percentAcc).*(sign(accelCostTemp  - tuning.max_accel*tuning.percentAcc)+1));                       # minus the max_vel and square
    #Add a weight
    accelCost *= tuning.accelWeight;

    #Create the soft costs? assuming soft costs are at the end
    softCosts = 0;
    #Make sure that there are enough weights
    tuning.softConstrWeights = checkQcoeffs(tuning.softConstrWeights, size(prob.PconstraintSoft,2))
    for i = 1:dim
        for j = 0:size(prob.PconstraintSoft,2)-1
            softCosts += (d[i,end-j]-prob.PconstraintSoft[i,end-j])^2 * tuning.softConstrWeights[end-j];
        end
    end
    
    #Create costmap costs
    #Note that the mapCost with a cell travel cost of 50 is on the order of 1000s
    mapCost = tuning.obstacleWeight*sum(prob.costmap[sol.cells]);

    #Return the total cost
    #println("Deriv: $derivCost, VelcCost: $velCost SoftCost: $softCosts Mapcost: $mapCost AccelCost: $accelCost")
    cost = derivCost + velCost + softCosts + mapCost+accelCost;
    return cost[1];
end

## form_df
#Description - forms the free constraints in a form that the cost function needs for good optimization
#Input
# prob - a PathProblem object with: 
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
#Output
# reshapedFree - the free constraints in the vector order of x's, y's, and z's
function form_df(prob)
    #Initialize to the expected size
    reshapedFree = zeros((2+prob.isDim3)*size(prob.PconstrFree,2))
    for i = 1:size(prob.PconstrFree,1)
        #Insert the correct values
        reshapedFree[(1:size(prob.PconstrFree,2))+(i-1)*size(prob.PconstrFree,2)] = prob.PconstrFree[i,:];
    end
    return reshapedFree;
end

## decompose_df
#TODO - make a void with reference parameters
#Description - the inverse of form_df; it forms the free constraints in the form of a problem's PconstrFree
#Input
# dF - the form of the free constraints returned by form_df
# prob - the PathProblem with:
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#Output
# shapedFree - the problem's PconstrFree variable with the optimized values
function decompose_df(dF, prob)
    #Create a dim variable
    dim = size(prob.PconstrFree,1)
    #Initialize the free constraint  variable
    shapedFree = zeros(dim, size(prob.PconstrFree,2));
    for i = 1:dim
        #Insert the correct values opposite the way in form_df
        shapedFree[i,:] = dF[(1:size(prob.PconstrFree,2))+(i-1)*size(prob.PconstrFree,2)];
    end 
    return shapedFree;
end

## createCostMap
#TODO - make it so objects are not at the beginning or end
#Descritption - creates a cost map for a problem
#Inputs
# obstacles - the number of random obstacles to make
# prob - the PathProblem object with
#   start_config::Array{Float64,2}   # Initial conditions with pos, vel, acc, jerk, ...
#   end_config::Array{Float64,2}     # Final conditions with pos, vel, acc, jerk,... in that order
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
#   grid_extent::Float64             # The extent of the cost map in meters
#   grid_resolution::Float64         # The resolution of the cost map in meters
#Output
# costmap - a costmap with obstacles that is 3D with one page only for 2D at the moment
function createCostMap(obstacles::Int64,prob::PathProblem)
    #Numbers to set
    costOfEmptyCell = 0;
    #Create an element number value
    n = round(Int64,ceil(prob.grid_extent/prob.grid_resolution));
    #Initialize to max size
    costmap = zeros(n,n,n)+costOfEmptyCell;
    #Create Objects
    for i = 1:obstacles
        #Random object centers 
        index1 = round(Int64,ceil(3/prob.grid_resolution))#round(rand()*n);
        index2 = round(Int64,ceil(3/prob.grid_resolution))#round(rand()*n);
        #Loop through and create
        for l = 1:size(costmap,1)
            for p = 1:size(costmap,2)
                costmap[l,p,1] += 900/(sqrt((l-index1)^2+(p-index2)^2));
                if(costmap[l,p,1]>255)
                    costmap[l,p,1] = 255;
                end
            end
        end
    end
    return costmap;
end

## debugPlotDash
#Description - plots the path and various derivative magnitudes as dashes with no background and can be used even if out of bounds
#Inputs
# sol - PathSol object with:
#   totTime::Float64            # Total time of the polynomial path in seconds
#   coeffs::Array{Float64,2}    # X, Y, Z, and Yaw coefficients
# prob - PathProblem object with:
#   isDim3::Bool                # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
# tuning - TuningParams object with:
#   timeRes::Int64              # How many points for determining limits and constraints
#Outputs
# Plots... with labels son! Aight! Cool...
function debugPlotDash(sol, prob, tuning)
    #Create times according to resolution
    plotTimes = linspace(0,sol.totTime,tuning.timeRes)
    #First figure is position
    figure(1)
    if(prob.isDim3)
        plot3D(   evaluate_poly(sol.coeffs[1,:],0,plotTimes), #x
                evaluate_poly(sol.coeffs[2,:],0,plotTimes), #y
                evaluate_poly(sol.coeffs[3,:],0,plotTimes), #z
                linestyle = "--");
        #Extra label for 3D: rad!
        zlabel("Z (m)")
    else
        plot(   evaluate_poly(sol.coeffs[1,:],0,plotTimes), #x
                evaluate_poly(sol.coeffs[2,:],0,plotTimes), #y 
                linestyle = "--");
    end
    #Labels Homeboi!
    title("Position")
    xlabel("X (m)")
    ylabel("Y (m)")
    #Titles and ylabels for reference in the for loop boi!
    titles = ["Velocity Magnitude", "Acceleration Magnitude", "Jerk Magnitude"]
    ylabels = ["Velocity Magnitude (m/s)", "Acceleration Magnitude(m/s^2)", "Jerk Magnitude (m/s^3)"]
    for i =1:3
        figure(i+1)
        if(prob.isDim3)
            #Various derivatives are sqaured summed and sqrt'ed
            plot(   plotTimes, sqrt(evaluate_poly(sol.coeffs[1,:],i,plotTimes).^2 + #x
                                    evaluate_poly(sol.coeffs[2,:],i,plotTimes).^2 + #y
                                    evaluate_poly(sol.coeffs[3,:],i,plotTimes).^2), #z
                                    linestyle = "--");
        else
            plot(   plotTimes, sqrt(evaluate_poly(sol.coeffs[1,:],i,plotTimes).^2 + #x
                                    evaluate_poly(sol.coeffs[2,:],i,plotTimes).^2), #y
                                    linestyle = "--");
        end
        title(titles[i])
        ylabel(ylabels[i])
        #Time label is same for all plots
        xlabel("Time (s)")
    end
end

## createRandomRestart
#Description - creates random values for the free constraints to start at before and optimization after a path has run 
# into an obstacle on the first pass
#Assumptions
# Max changes in derivatives are provided in tuning
#Inputs
# prob - a PathProblem with:
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   PconstraintOrders::Array{Int64,1}# The order corresponding to each constraint in the total constriant vector
# tuning - a TuningParams object with:
#   max_vel::Float64                   # Max total velocity the path is restricted to
#   max_accel::Float64                 # Max total acceleration the path is restricted to
#   max_jerk::Float64                  # Max total jerk the path is restricted to
#Outputs
# randConstrFree - an updated PconstrFree
function createRandomRestart(prob, tuning)
    #Create a dimension variable
    dim = 2+ prob.isDim3;
    #Initialize 
    randConstrFree = prob.PconstrFree;
    #Loop throught the row of randConstrFree
    for i = 1:size(randConstrFree,1);
        #Loop through each constraint order
        for j = 1:size(prob.PconstrFree,2)
            #if statements for what value to put in
            if(prob.PconstraintOrders[j+size(prob.PconstrFixed,2)] == 0)
                #For position points
                #max value to add, times 2, make rand() -0.5 to 0.5 and add the position that is desired
                randConstrFree[i,j] = tuning.posMaxAdd*2*(0.5-rand())+prob.PconstrFree[i,j]
            elseif(prob.PconstraintOrders[j+size(prob.PconstrFixed,2)]==1)
                #For vel points
                #max value to add, times 2, divide by dimensions, make rand() -0.5 to 0.5
                randConstrFree[i,j] = tuning.velMaxAdd*2*(0.5-rand())+prob.PconstrFree[i,j]#tuning.max_vel*2/dim*(0.5-rand());
            elseif(prob.PconstraintOrders[j+size(prob.PconstrFixed,2)]==2)
                #For accel points
                #max value to add, times 2,divide by dimensions, make rand() -0.5 to 0.5
                randConstrFree[i,j] = tuning.accelMaxAdd*2*(0.5-rand())+prob.PconstrFree[i,j]#tuning.max_accel*2/dim*(0.5-rand());
            else
                #For everything else
                #max value to add, times 2,divide by dimensions, make rand() -0.5 to 0.5
                randConstrFree[i,j] = tuning.jerkMaxAdd*2*(0.5-rand())+prob.PconstrFree[i,j]#tuning.max_jerk*2/dim*(0.5-rand());
            end
        end
    end

    return randConstrFree;
end

## pathCellPlot
#Description - for the solution, the path is plotted and the cells it goes through set to 255 to be highlighted on a
# plot with the costmap
#Assumption
# Cells are within the costmap
# The path is final
# Costmap is 2D
#Inputs
# prob - a PathProblem with:
#   costmap::Array{Float64,3}        # The 3D voxel occupancy grid
#   grid_extent::Float64             # The extent of the cost map in meters
# sol - a PathSol with:
#   totTime::Float64            # Total time of the polynomial path in seconds
#   coeffs::Array{Float64,2}    # X, Y, Z, and Yaw coefficients
#   cells::Array{Int64,1}       # The indices of the costmap that the path goes through
# tuning - TuningParams with
#   timeRes::Int64                     # How many points for determining limits and constraints
#Outputs
# Plot of the cost map, the path, and the highlighted cells
function pathCellPlot(prob, sol, tuning)
    #Create times according to resolution
    plotTimes = linspace(0,sol.totTime,tuning.timeRes)
    #Read in costmap
    costmap = prob.costmap;
    #Create this figure in figure 5
    figure(5)
    subplot(1,2,1)
    imshow( flipCostmap(costmap),
            cmap = "gray", 
            interpolation="none",
            extent=[0,prob.grid_extent,0,prob.grid_extent])
    if(prob.isDim3)
        plot3D(   evaluate_poly(sol.coeffs[1,:],0,plotTimes), #x
                evaluate_poly(sol.coeffs[2,:],0,plotTimes), #y
                evaluate_poly(sol.coeffs[3,:],0,plotTimes));#z
        #Extra label for 3D: rad!
        zlabel("Z (m)")
    else
        plot(   evaluate_poly(sol.coeffs[1,:],0,plotTimes), #x
                evaluate_poly(sol.coeffs[2,:],0,plotTimes));#y 
    end
    #Labels Homeboi!
    title("Final Path")
    xlabel("X (m)")
    ylabel("Y (m)")
    subplot(1,2,2)
    #Add the cells to the costmap! whooo!
    for i in sol.cells
        costmap[i]=255;
    end
    #Repeat the above
    imshow( flipCostmap(costmap),
            cmap = "gray", 
            interpolation="none",
            extent=[0,prob.grid_extent,0,prob.grid_extent])
    if(prob.isDim3)
        plot3D(   evaluate_poly(sol.coeffs[1,:],0,plotTimes), #x
                evaluate_poly(sol.coeffs[2,:],0,plotTimes), #y
                evaluate_poly(sol.coeffs[3,:],0,plotTimes));#z
        #Extra label for 3D: rad!
        zlabel("Z (m)")
    else
        plot(   evaluate_poly(sol.coeffs[1,:],0,plotTimes), #x
                evaluate_poly(sol.coeffs[2,:],0,plotTimes));#y 
    end
    #Labels Homeboi!
    title("Final Path with Cells")
    xlabel("X (m)")
    ylabel("Y (m)")
    subplot(1,2,2)
end

## selfGradientDescent
#Description - conduct a gradient descent optimization over the free constraints of a PathProblem for a 
# set amount of iterations
#Assumptions
# All constraints are perturbed the same amount
#Input
# sol - a PathSol with:
#   totTime::Float64            # Total time of the polynomial path in seconds
#   coeffs::Array{Float64,2}    # X, Y, Z, and Yaw coefficients
#   cells::Array{Int64,1}       # The indices of the costmap that the path goes through
#   cost::Float64               # The total cost of the path
# prob - the PathProblem with the free constraints to be:
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
# tuning - the TuningParams object with the following for this portion:
#   precision::Float64                 # The tolerance to solve within
#   iterations::Int64                  # Number of times to do optimization step
#   perturbation::Float64              # How much to start out changing free variables by
# solvOpt - a PolyPathSolver object with everything
#Output
# soln - the optimized PathSol object
# PconstrFree - the optimizing free constraints
# outOfBounds - a flag of whether the path went out of bounds or not
function selfGradientDescent(sol::PathSol,prob::PathProblem,tuning::TuningParams,solvOpt::PolyPathSolver)
    #Start optimization as unoptimized
    #Create a copy of the problem
    oldCost = sol.cost;
    soln = PathSol(sol.totTime, sol.coeffs, sol.cells, sol.cost);
    probOpt = prob;
    perturb = tuning.perturbation;
    outOfBounds = false;
    rowFree = size(prob.PconstrFree,1);
    colFree = size(prob.PconstrFree,2);
    counter = 0;
    rateChange = zeros(rowFree*colFree);
    while(counter<tuning.iterations)
        #Start while?
        #For each free variable perturb by an amount
        for i = 1:rowFree
            for j = 1:colFree
                #Hold the pre perturbed stuff
                #Perturb specific free variable the amount
                probOpt.PconstrFree[i,j] += perturb;

                #solve for the polynomial coefficients
                soln.coeffs = solvePolysInitially(probOpt,solvOpt);

                #Check in bounds and collect the cells
                soln.cells, outOfBounds = occupancyCellChecker(soln, probOpt, tuning);

                #Break if out of bounds and display an error
                if(outOfBounds)
                    println("Plan Fail: Went out of Bounds in First Optimize Step")
                    #Add a return here
                    #return (PathSol(soln.totTime, soln.coeffs, soln.cells, soln.cost), probOpt.PconstrFree, outOfBounds);
                end

                #Create a holder for the free constr
                dF = form_df(probOpt);

                #Calculate the costs
                costNew = costFunc(dF, soln, probOpt, solvOpt, tuning)
                
                #Record the change in cost
                rateChange[(i-1)*colFree+j] = (costNew-oldCost)/perturb;

                #Reset the free constraint so it doesn't affect other ones
                probOpt.PconstrFree[i,j] -= perturb
            end    
        end
        
        if(all(rateChange .== 0.0))
            gradDir = zeros(size(rateChange));
        else
            #Only Normalize the rate of change vector if possible
            gradDir = normalize!(rateChange)
        end
        println(gradDir)

        #Perturb in the negative gradient direction
        for i = 1:rowFree
            for j = 1:colFree
                #Perturb in the negative gradient direction
                probOpt.PconstrFree[i,j] += perturb*-gradDir[(i-1)*colFree+j];
            end    
        end

        #solve for the polynomial coefficients
        soln.coeffs = solvePolysInitially(probOpt,solvOpt);
        
        #Check in bounds and collect the cells
        soln.cells, outOfBounds = occupancyCellChecker(soln, probOpt, tuning);

        #Break if out of bounds and display an error
        if(outOfBounds)
            println("Plan Fail: Went out of Bounds In Optimization")
            #Add a return here
            #return (PathSol(soln.totTime, soln.coeffs, soln.cells, soln.cost), probOpt.PconstrFree, outOfBounds);
        end

        #Create a holder for the free constr
        dF = form_df(probOpt);

        #Calculate the costs
        costNew = costFunc(dF, soln, probOpt, solvOpt, tuning)

        #Check if step is too small or if change is too small and break if within precision
        if(all(abs(gradDir*perturb) .< tuning.precision) || abs(costNew - oldCost) < tuning.precision)
            #Exit optimization with a return
            return (PathSol(soln.totTime, soln.coeffs, soln.cells, soln.cost), probOpt.PconstrFree, outOfBounds);
        elseif(costNew - oldCost < 0)
            #If step was a decrease step keep that and increase step size
            perturb *= 2;
            #Update cost
            oldCost = costNew;
            soln.cost = oldCost;
            
        else
            #If step was an increase revert the step and half the step size
            for i = 1:rowFree
                for j = 1:colFree
                    #Reverse perturb in the negative gradient direction
                    probOpt.PconstrFree[i,j] -= perturb*-gradDir[(i-1)*colFree+j];
                end    
            end
            perturb /= 2;
        end
        #Update counter
        counter += 1;
    end

    #If made it here went past all the iterations
    println("Optimization went through all the iterations")
    
    #Return the solution the updated free constraints and the out of bounds flag
    return (PathSol(soln.totTime, soln.coeffs, soln.cells, soln.cost), probOpt.PconstrFree, outOfBounds);
end
