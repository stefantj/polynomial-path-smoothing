## constructConstr!
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
    constrFixed =  zeros(2+prob.isDim3,degree-length(softIndex)-prob.dof);
    for i = 1:(2+prob.isDim3)
        constrFixed[i,:] = [prob.start_config[i,:]' prob.end_config[i,notSoftIndex]'];
    end
    #Create the free vector with soft constraints at the end
    constrFree =  zeros(2+prob.isDim3,length(softIndex)+prob.dof);
    for i = 1:(2+prob.isDim3)
        constrFree[i,:] = [zeros(size(1:prob.dof))' constrSoft[i,:]']
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
    #Maxsure that if max_vel_vec is zero you do not get NAN
    if(max_vel_vec != zeros(2+prob.isDim3))
        max_vel_vec = normalize!(max_vel_vec)*max_vel
    end
    #Add as a soft constraint being careful about what exists
    if(size(prob.end_config,2) <= 1)
        prob.end_config = [prob.end_config[1:(2+prob.isDim3),1] max_vel_vec];
        prob.soft_constr = [prob.soft_constr[1]; true];
    else
        prob.end_config[1:(2+prob.isDim3),2] = max_vel_vec;
        prob.soft_constr[2] = true;
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

## udpateFreeConstr
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
        prob.PconstrFree[i,:] = solvStuff.PoptimizeMatrix * prob.PconstrFree[i,:] + addSofts[i,:];
    end
    return prob.PconstrFree;
end


## solvePolysInitially
#Description - will do the matrix multiplication of optimizeMat and the fixed constraints to get the optimal free constraints
# w.r.t. the derivative costs.
#Assumptions
# All the dimensions of the problem have the same number of soft constraints and free constraints
#Inputs
# prob - the PathProblem object with the following values:
#   PconstrFixed::Array{Float64,2}   # All fixed constraints
#   PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
#   isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3  
#   Pdegree::Int64                   # The degree of the polynomial to be made
# solvStuff - a PolyPathSolver object with the following values:
#   PA_inv::Array{Float64,2}         # The matrix that transforms constraints to coefficients A_inv*d=p
#Outputs
# coeffs - the polynomial coefficients in a matrix
function solvePolysInitially(prob,solvStuff)
    #initialize the coeeffs variable
    coeffs = zeros(2+prob.isDim3,prob.Pdegree)
    #for each dimension
    for i = 1:(2+prob.isDim3)
        #p = A^-1 * d
        coeffs[i,:] = solvStuff.PA_inv * [prob.PconstrFixed[i,:]' prob.PconstrFree[i,:]']';
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