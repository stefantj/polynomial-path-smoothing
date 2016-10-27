using PyPlot # So plots can be used
using JuMP   # Setup using JuMP so that we have access to its methods
using Optim  # For optim.jl use

#Classes/Objects without the functions
#Convention Notes: 
# In the future try to use matrices instead of separate arrays for similar values
# For matrices: x,y,z,yaw -> 1,2,3,4
#               pos,vel,accel,jerk,...->1,2,3,4,...
# Using matrices would require all polynomials to be the same degree and exact matching of values 
# A capital P will mean that a variable is private and should not be worried for initialization

#ASSUMPTION: all paths start at t = 0 seconds
type PathSol                    # The polynomial path solution
    totTime::Float64            # Total time of the polynomial path in seconds
    coeffs::Array{Float64,2}    # X, Y, Z, and Yaw coefficients
    cells::Array{Int64,1}       # The indices of the costmap that the path goes through
    cost::Float64               # The total cost of the path
end

#ASSUMPTION: the path problem assumes each dimension is solved in the same way
type PathProblem                     # The parameters that define the problem
    start_config::Array{Float64,2}   # Initial conditions with pos, vel, acc, jerk, ...
    end_config::Array{Float64,2}     # Final conditions with pos, vel, acc, jerk,... in that order
    soft_constr::Array{Bool,1}       # A vector that determines whether the end configurations are soft or not
    DijkstraNotFMT::Bool             # A flag to determine what how to solve the problem
    PconstrFixed::Array{Float64,2}   # All fixed constraints
    PconstrFree::Array{Float64,2}    # A holder for all soft and free contraints
    PtimeIndex::Array{Int64,1}       # The vector of the times at which the constraints apply
    PconstraintOrders::Array{Int64,1}# The order corresponding to each constraint in the total constriant vector
    PconstraintSoft::Array{Float64,2}# Soft constraints ordered by ascending order
    isDim3::Bool                     # The flag for 3D; 2 will be added later so the dimension is only 2 or 3
    dof::Int64                       # The extra degrees that would be completely free to optimize over
    costmap::Array{Float64,3}        # The 3D voxel occupancy grid
    grid_extent::Float64             # The extent of the cost map in meters
    grid_resolution::Float64         # The resolution of the cost map in meters
    Pgrid_elementNum::Int64          # The max index of the grid of the cost map
    Pdegree::Int64                   # The degree of the polynomial to be made
end

type TuningParams                      # All the knobs to turn for tuning
    q_coeff::Array{Int64,1}            # The weights for the derivatives to minimize with the Q matrix
    softConstrWeights::Array{Float64,1}# The vector of weights to make soft constraints harder
    numericVelWeight::Float64          # The weight for the cost for not being at max_vel at any time
    obstacleWeight::Float64            # The weight for the cost of going through a cell in the costmap
    derivativeWeight::Float64          # The weight for the combined derivative cost term
    precision::Float64                 # The tolerance to solve within
    iterations::Int64                  # Number of times to do optimization step
    timeStep::Float64                  # The increment in time to determine the cells that a path goes through
    timeStart::Float64                 # Time to start optimization at
    aggressParam::Float64              # Closer to 0 means check every point on the path, infinity => check no points
    timeIncrease::Float64              # The increment by which to increase time if path is found infeasible
    max_vel::Float64                   # Max total velocity the path is restricted to
    max_accel::Float64                 # Max total acceleration the path is restricted to
    max_jerk::Float64                  # Max total jerk the path is restricted to
    numberOfRandomRestarts::Int64      # Number of random restarts to do before giving up
    timeWeight::Float64                # Weight applied to time
    timeRes::Int64                     # How many points for determining limits and constraints
    precisionVel::Float64              # Precision in all derivatives that it is okay to be over limits by
    percentAcc::Float64                # The decimal percent by which to bias the path to accelerate at
    accelWeight::Float64               # Weight on acceleration term to bias path to accelerate near that
end


type PolyPathSolver                    # Important things for solving for a polynomial path
    PA_inv::Array{Float64,2}           # The matrix that transforms constraints to coefficients A_inv*d=p
    PQ::Array{Float64,2}               # Matrix representation of the q_coeff weights
    PoptimizeMatrix::Array{Float64,2}  # Gives the optimal free/soft constraints given fixed constraints
    PoutOfBounds::Bool                 # True if path ever goes outside the costmap scope
    PunVerified::Bool                  # True if path is not feasible based on derivative restrictions
    counterRestart::Int64              # Counts how many restarts have been done
    counterVerified::Int64             # Counts how many times the time has been increased to make path feasible
end




##Position Lotus Test
#Description - Start and end points are separated an initial distance from each other. The path is tested for
# this separation in directions that are determined from a fraction of a 180 rotation. For two more times, the 
# distance is increased and the test run again. The start and end are centered around (0,0).
#Inputs
# initial_distance - the initial distance between the start and end points
# tuningParams - the tuning parameter configuration you want to test
# number2Rotate  - the number of times that the positions are rotated for a certain distance
# distanceChange - the amount distance should change by
# figureNum - the figure number on which the plot will print
# pathProb - the defining variables of the path planning problem
#Expected Outcome - A graph of straight paths that go through (0,0) with dots at the start and end points. 
# It will look like a flower.
function lotusTest(initial_distance::Float64, 
                   tuningVars, #Will be of the tuningParam type
                   number2Rotate::Int64,
                   distanceChange::Float64, 
                   figureNum::Int64,
                   pathProb) #::PathProblem)
    #Record the initial_distance divided by two in a value that will increase later
    distanceHolder = initial_distance/2;

    #Create a vector of angles to test at
    angles2Test = collect(0:(pi/(number2Rotate+1)):pi)
    #Create figure
    figure(figureNum)
    #For 3 times
    for numberOFTests = 1:3
        #Create the matrix of start and end points x = 1st column and y = 2nd column
        startPoints = [-distanceHolder*cos(angles2Test) -distanceHolder*sin(angles2Test)]
        endPoints = [distanceHolder*cos(angles2Test) distanceHolder*sin(angles2Test)]
        #For each set of points
        for looper = 1:size(startPoints, 1)
            #########Call the planner###################################
            #TODO make this more efficient
            #Create a copy of the pathProb and 
            probCopy = pathProb;
            #This is a 2D test at the moment
            probCopy.isDim3 = false;
            #make sure the input is of the right size
            if(size(probCopy.start_config,1) > 2)
                probCopy.start_config = probCopy.start_config[1:2,:];
            end
            if(size(probCopy.end_config,1) > 2)
                probCopy.end_config = probCopy.end_config[1:2,:];
            end
            #change the initial and final configurations
            probCopy.start_config[:,1] = [startPoints[looper,1]; #x start position
                                          startPoints[looper,2]];#y
            
            probCopy.end_config[:,1] = [endPoints[looper,1]; #x end position
                                        endPoints[looper,2]];#y
            #Call the planner
            solution = runPathPlanner(tuningVars,probCopy);
            #plot each iteration
            lotusTimes = linspace(0,solution.totTime,100);
            plot(evaluate_poly(solution.coeffs[1,:],0,lotusTimes),#x
                 evaluate_poly(solution.coeffs[2,:],0,lotusTimes))#y
            ############################################################
            #scatter plot the endpoints so they show up
            scatter(startPoints[looper,1], #x
                    startPoints[looper,2], #y 
                    color = "red") 
            scatter(endPoints[looper,1], #x
                    endPoints[looper,2], #y
                    color = "green") 
        end
        #Increase the distance by half the distanceChange since half is added to each side
        distanceHolder += distanceChange/2
    end

    #title for entire thing
    title("Lotus Test: Starts are Red and Ends are Green")
end

##Flip Costmap for Plotting Function
#Description: Flips a matrix of a costmap so that it can be shown correctly on imshow
#Inputs
# costmap2Flip - a 2D costmap to flip
#Outputs
# flippedCostmap - the flipped costmap
function flipCostmap(costmap)
    #Create a the flipping costmap
    flippedCostmap = zeros(size(costmap,1),size(costmap,2))
    for h =1:size(costmap,1)
        for g = 1:size(costmap,2)
            #Flip the columns
            flippedCostmap[h,g] = costmap[h,size(costmap,2)-g+1];
        end
    end
    #Switch rows and columns with transpose 
    flippedCostmap=flippedCostmap';
    #Return the costmap
    return flippedCostmap;
end

##Easy Object in the Way Function
#Description - start and end points are chosen at (1,5) and (9,5) with a circular object centerd around (5,5).
# The planner must then plan around the object without hitting it. Path will then be assessed common sensically.
#Inputs
# bumpTuning - the tuning parameters for the planner
# pathProblem - the characteristics of the planning problem
# bumpSize - a nonegative integer that is less than or equal to 5000
# figNum - the number of the figure to print to
# bumpYLoc - the choice of where the bump is located along the axis x = 5
#Expected Success Outcomes - The path bends up or down to avoid the obstacle as if stretched by the obstacle
# without hitting it
#Expected Failure Modes
# Running into the obstacle - not enough optimization away from obstacle
# Going out of bounds - not understood why but this is stupid so work out
# Overavoiding the obstacle - too much optimization away from obstacle
function easyAvoidBumpTest(bumpTuning::TuningParams, pathProblem::PathProblem, bumpSize::Int64, figNum::Int64 = 1,
    bumpYLoc::Float64 = 5.0)
    #Make sure bump size is acceptable
    if(bumpSize > 5000) #5000 was just a feely number nothing really behind it
        println("Test Abort: $bumpSize is greater than 5000")
        return -1;
    end
    #TODO: Find out how not to do copying safely and reliably
    #Copy the problem locally
    bumpProb = pathProblem;
    #Set the start
    bumpProb.start_config[1,1] = 1; #x pos
    bumpProb.start_config[2,1] = 5; #y pos
    #Set the end
    bumpProb.end_config[1,1] = 9; #x pos
    bumpProb.end_config[2,1] = 5; #y pos
    #Create a costmap that has a bump in the middle
    #Create the length of the map and update problem accordingly
    bumpProb.grid_extent = 10.0;
    bumpProb.grid_resolution = 0.1;
    bumpProb.Pgrid_elementNum = round(Int64,ceil(bumpProb.grid_extent/bumpProb.grid_resolution));
    bumpMap = zeros(bumpProb.Pgrid_elementNum,
                    bumpProb.Pgrid_elementNum,
                    bumpProb.Pgrid_elementNum);
    #TODO: Make the order so that if there is an error some calcs are not done
    #TODO: Make it so that any grid extent can be used
    #Make sure bump y location is withing the extents of the grid
    if(bumpYLoc > bumpProb.grid_extent || bumpYLoc < 0.0)
        println("Test Abort: $bumpYLoc is out of the bounds [0.0,$(bumpProb.grid_extent)]")
        return -1;
    end
    #At the center add a bump
    centerX = round(bumpProb.Pgrid_elementNum/2);
    centerY = round(Int64,bumpYLoc/bumpProb.grid_resolution);
    for row = 1:size(bumpMap,1)
        for col = 1:size(bumpMap,2)
            #Lazy way of creating a bump TODO: Make not lazy
            bumpMap[row,col,1] += bumpSize/(sqrt((row-centerX)^2+(col-centerY)^2)+1);
            #Round large values to 255
            if(bumpMap[row,col,1]>255)
                bumpMap[row,col,1] = 255;
            end
        end
    end
    #########Call the planner###################################
    #Call the planner
    solution = runPathPlanner(bumpTuning,bumpProb);
    #plot each iteration
    bumpTimes = linspace(0,solution.totTime,100);
    figure(figNum)
    plot(evaluate_poly(solution.coeffs[1,:],0,bumpTimes),#x
         evaluate_poly(solution.coeffs[2,:],0,bumpTimes),#y
         color = "yellow")
   ############################################################
    #Check if bump is there
    imshow( flipCostmap(bumpMap),
            cmap = "gray", 
            interpolation="none",
            extent=[0,bumpProb.grid_extent,0,bumpProb.grid_extent])
    #Don't Forget the title
    title("Easy to Avoid Obstacle Test")
end

##Doorway Test 
#Description - start and end points are chosen at (1,5) and (9,5) with a doorway centerd on x = 5 and a user
# specified height location and width. The planner must then plan around the object without hitting it. 
# Path will then be assessed common sensically.
#Inputs
# wallTuning - the tuning parameters for the planner
# pathProblem - the characteristics of the planning problem
# wallSize - a nonegative integer that is less than or equal to 5000
# figNum - the number of the figure to print to
# doorYLoc - the choice of where the door is located along the axis x = 5
# doorWidth - the width of the door, must be less than the grid extent
# wallBehindDoor - a flag that will include a wall behind
#Expected Success Outcomes - The path bends up or down to go in the doorway's middle
#Expected Failure Modes
# Running into the obstacle - not enough optimization away from obstacle
# Going out of bounds - not understood why but this is stupid so work out
function doorwayTest(wallTuning::TuningParams, pathProblem::PathProblem, wallSize::Int64, figNum::Int64, 
    doorYLoc::Float64, doorWidth::Float64, wallBehindDoor::Bool = false, wallDistance::Float64 = 2.0)
    #Check if the door width is positive.
    if(doorWidth <= 0)
        println("Test Abort: Invalid door width; it must be positive")
        return -1;
    end
    #TODO: try not to copy
    wallProb = problem;
    #Set the start
    #TODO: Move start and ends to user specified values
    wallProb.start_config[1,1] = 1; #x pos
    wallProb.start_config[2,1] = 5; #y pos
    #Set the end
    wallProb.end_config[1,1] = 9;
    wallProb.end_config[2,1] = 5;
    #Create the length of the map and update problem accordingly
    wallProb.grid_extent = 10.0;
    wallProb.grid_resolution = 0.1;
    wallProb.Pgrid_elementNum = round(Int64,ceil(wallProb.grid_extent/wallProb.grid_resolution));
    #Create a costmap that has a wall in the middle
    wallMap = zeros(wallProb.Pgrid_elementNum,
                    wallProb.Pgrid_elementNum,
                    wallProb.Pgrid_elementNum);
    #Create the start and end of the door indeces
    doorStart = round(Int64,doorYLoc/wallProb.grid_resolution)-
                round(Int64,doorWidth/2);
    doorEnd = round(Int64,doorYLoc/wallProb.grid_resolution)+
              round(Int64,doorWidth/2);
    #Make sure the door doesn't over index
    if(doorStart < 1)
        doorStart = 1;
    end
    if(doorEnd > wallProb.Pgrid_elementNum)
        doorEnd = wallProb.Pgrid_elementNum;
    end
    #TODO: Come up with a better way to create a wall in the costmap
    #At the center add a doorway by adding successive bumps at locations along x = 5 that are not in the door
    for wallLoop in [1:doorStart; doorEnd:size(wallMap,1)]
        centerX = round(Int64, wallProb.Pgrid_elementNum/2);
        centerY = wallLoop;
        for row = 1:size(wallMap,1)
            for col = 1:size(wallMap,2)
                #Lazy way of creating a wall TODO: Make not lazy
                wallMap[row,col,1] += wallSize/(sqrt((row-centerX)^2+(col-centerY)^2));
                #Round large values to 255
                if(wallMap[row,col,1]>255)
                    wallMap[row,col,1] = 255;
                end
            end
        end
    end
    #Add the wall if desired
    if(wallBehindDoor)
        #TODO: Come up with a better way to create a wall in the costmap
        #At the center add a doorway by adding successive bumps at locations along x = 5 that are not in the door
        for wallLoop = doorStart:doorEnd
            centerX = round(Int64,  wallProb.Pgrid_elementNum/2+
                                    wallDistance/wallProb.grid_resolution);
            centerY = wallLoop;
            for row = 1:size(wallMap,1)
                for col = 1:size(wallMap,2)
                    #Lazy way of creating a wall TODO: Make not lazy
                    wallMap[row,col,1] += wallSize/(sqrt((row-centerX)^2+(col-centerY)^2));
                    #Round large values to 255
                    if(wallMap[row,col,1]>255)
                        wallMap[row,col,1] = 255;
                    end
                end
            end
        end
    end
    #########Call the planner###################################
    #Call the planner
    solution = runPathPlanner(wallTuning,wallProb);
    #plot each iteration
    wallTimes = linspace(0,solution.totTime,100);
    figure(figNum)
    plot(evaluate_poly(solution.coeffs[1,:],0,wallTimes),#x
         evaluate_poly(solution.coeffs[2,:],0,wallTimes),#y
         color = "yellow")
    ############################################################
    #Check if wall is there
    imshow( flipCostmap(wallMap),
            cmap = "gray", 
            interpolation="none",
            extent=[0,wallProb.grid_extent,0,wallProb.grid_extent])
    #Create a title
    if(wallBehindDoor)
        title("Psuedo Doorway Test: With Wall")
    else
        title("Psuedo Doorway Test: Without Wall")
    end
end

##Hallway test
#Description - start and end points are chosen at (1,4) and (9,6) with  walls at y = 2 and 8. The planner is then
# invoked to create an optimized path that will be judged common sensically.
#Inputs
# hallTuning - the tuning parameters for the planner
# pathProblem - the characteristics of the planning problem
# wallSize - a nonegative integer that is less than or equal to 5000
# figNum - the number of the figure to print to
#Expected Success Outcomes - The path is close to a straight line that can bend towards the middle of the hallway
#Expected Failure Modes
# Unoptimal path - consult common sensitist for direction
# Going into wall - not understood why but this is stupid so work out
function hallwayTest(hallTuning::TuningParams, pathProblem::PathProblem, wallSize::Int64, figNum::Int64)
    hallProb = problem;
    #TODO: get rid of the hardcoded numbers for user input
    centerY2 = round(Int64,2/hallProb.grid_resolution);
    centerY8 = round(Int64,8/hallProb.grid_resolution);
    #Set the start
    #TODO: Move start and ends to user specified values
    hallProb.start_config[1,1] = 1;
    hallProb.start_config[2,1] = 4;
    #Set the end
    hallProb.end_config[1,1] = 9;
    hallProb.end_config[2,1] = 6;
    #Create the length of the map and update problem accordingly
    hallProb.grid_extent = 10.0;
    hallProb.grid_resolution = 0.1;
    hallProb.Pgrid_elementNum = round(Int64,ceil(hallProb.grid_extent/hallProb.grid_resolution));
    #Create a costmap that has a wall in the middle
    hallMap = zeros(hallProb.Pgrid_elementNum,
                    hallProb.Pgrid_elementNum,
                    hallProb.Pgrid_elementNum);
    #TODO: Come up with a better way to create a wall in the costmap
    #Add two walls at y = 2 and 8 
    for hallLoop = 1:size(hallMap,1)
        centerX = hallLoop;
        for row = 1:size(hallMap,1)
            for col = 1:size(hallMap,2)
                #Lazy way of creating a wall TODO: Make not lazy
                #The wall at y=2
                hallMap[row,col,1] += wallSize/(sqrt((row-centerX)^2+(col-centerY2)^2));
                #Wall at y=8
                hallMap[row,col,1] += wallSize/(sqrt((row-centerX)^2+(col-centerY8)^2));
                #Round large values to 255
                if(hallMap[row,col,1]>255)
                    hallMap[row,col,1] = 255;
                end
            end
        end
    end
    #########Call the planner###################################
    #Call the planner
    solution = runPathPlanner(hallTuning,hallProb);
    #plot each iteration
    hallTimes = linspace(0,solution.totTime,100);
    figure(figNum)
    plot(evaluate_poly(solution.coeffs[1,:],0,hallTimes),#x
         evaluate_poly(solution.coeffs[2,:],0,hallTimes),#y
         color = "yellow")
    ############################################################
    #Check if wall is there
    imshow( flipCostmap(hallMap),
            cmap = "gray", 
            interpolation="none",
            extent=[0,hallProb.grid_extent,0,hallProb.grid_extent])
    #Create a title
    title("Hallway Test")
end


##Evaluate Poly Function
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


##Create a Function Stub for the path planner
function runPathPlanner(tuning::TuningParams, 
                        prob::PathProblem)
    seedForRandomness = 1;
    srand(seedForRandomness);

    #Initialize solver and solution
    solution = PathSol( 0.0,     # totTime::Float64 
                        [0.0]'', # coeffs::Array{Float64,2}  
                        [1],     # cells::Array{Int64,1}  
                        0.0)     # cost::Float64         
    solvHelp = PolyPathSolver(  [0.0]'', #PA_inv::Array{Float64,2}           
                                [0.0]'', #    PQ::Array{Float64,2}              
                                [0.0]'', #    PoptimizeMatrix::Array{Float64,2}  
                                true,    #    PoutOfBounds::Bool                 
                                true,    #    PunVerified::Bool                
                                10,      #    counterRestart::Int64             
                                100);    #    counterVerified::Int64    
    #Create a costmap
    prob.costmap = createCostMap(1,problem);

    #If dijkstras Create a normalized direction vector
    if(prob.DijkstraNotFMT)
        addDirectedSpeed!(prob, tuning.max_vel);
    end

    #Construct the constraint vectors, orders, and timeIndex given start, end, dof and soft constraint vectors
    #TODO: make void functions with pointers to avoid so much copying
    prob = constructConstr(prob)

    #initialize time to a reasonable time so polynomial doesn't explode, a reasonable time will be
    solution.totTime = initializeTime(prob, tuning);

    #while loop until verified 
    while(solvHelp.PunVerified)
        #Increment time on every loop
        solution.totTime += tuning.timeIncrease;

        #Form A_inv
        solvHelp.PA_inv = constr_Ainv(  prob.PconstraintOrders,
                                        prob.PtimeIndex*solution.totTime, #Multiply by the total time for proper function
                                        prob.Pdegree);

        #Add/subtract stuff to/from q_coeff if not equal to the degree
        tuning.q_coeff = checkQcoeffs(tuning.q_coeff, prob.Pdegree)

        #Form Q
        solvHelp.PQ = form_Q(tuning.q_coeff,solution.totTime);

        #Form OptimizeMat
        solvHelp.PoptimizeMatrix = form_OptimizeMat(prob,tuning,solvHelp);

        #update the free constraints as necessary
        prob.PconstrFree = updateFreeConstr(prob,solvHelp);

        #solve for the polynomial coefficients
        solution.coeffs = solvePolysInitially(prob,solvHelp);

        #Check in bounds and collect the cells
        solution.cells, solvHelp.PoutOfBounds = occupancyCellChecker(solution, prob, tuning);

        #Break if out of bounds and display an error
        if(solvHelp.PoutOfBounds)
            println("Plan Fail: Went out of Bounds")
            break;
        end

        #Verify good path
        errorVals, errorTypes = simpleVerifyFeas(solution, tuning)

        #If not verified increment time
        if(!isempty(errorVals))
            #Repeat loop over
            continue; 
        end
        #If verified do not say verified  until next verify
        
        #Create a holder for the free constr
        df = form_df(prob);

        #optimize free constraints with limits built in using the COST function
        #Do one iteration solve next iteration repeat
        for i = 1:tuningI.iterations
            result = optimize(x -> costFunc(x, solution, prob, solvHelp, tuning), df,
                                            GradientDescent(),
                                            OptimizationOptions(iterations = 1));

            #Save the result of the optimization
            df = Optim.minimizer(result)

            #Resolve df into the respective free constraints and solve again
            prob.PconstrFree = decompose_df(df, prob);

            #solve for the polynomial coefficients
            solution.coeffs = solvePolysInitially(prob,solvHelp);

            #Check in bounds and collect the cells
            solution.cells, solvHelp.PoutOfBounds = occupancyCellChecker(solution, prob, tuning);
        end


        #Break if out of bounds and display an error
        if(solvHelp.PoutOfBounds)
            println("Plan Fail: Went out of Bounds")
            break;
        end

        #Verify good path Number 2
        errorVals, errorTypes = simpleVerifyFeas(solution, tuning)

        #If not verified increment time
        if(!isempty(errorVals))
            #Repeat loop over
            continue; 
        end
        solvHelp.PunVerified = false;

        #end while loop
    end


    #Start random restarts if hitting an obstacle
    restartCounter = 0;
    while(any(prob.costmap[solution.cells] .>= 255) && restartCounter < tuning.numberOfRandomRestarts)

        #Print message about hitting an obstacle and random restarting
        println("Hit and Obstacle; Trying a random restart")

        #Create random start
        prob.PconstrFree = createRandomRestart(prob, tuning);

        #optimize - note there is no verification step anymore so no more increasing time this could cause problems
        # later
        #Create a holder for the free constr
        df = form_df(prob);

        #1)OPTIMIZE free constraints with limits built in; 2)COST function
        #Do one iteration solve next iteration repeat
        for i = 1:tuning.iterations
            result = optimize(x -> costFunc(x, solution, prob, solvHelp, tuning), df,
                                            GradientDescent(),
                                            OptimizationOptions(iterations = 1));

            #Save the result of the optimization
            df = Optim.minimizer(result)

            #Resolve df into the respective free constraints and solve again
            prob.PconstrFree = decompose_df(df, prob);

            #solve for the polynomial coefficients
            solution.coeffs = solvePolysInitially(prob,solvHelp);

            #Check in bounds and collect the cells
            solution.cells, solvHelp.PoutOfBounds = occupancyCellChecker(solution, prob, tuning);
        end

        #Increment the counter
        restartCounter += 1; 
    end

    #Print if failed to plan a path around obstacle
    if(restartCounter >= tuning.numberOfRandomRestarts)
        println("Path Planning Failed even with Restarts: Sad Face :(")
    end

    #Return path
    return solution;

    
end


##Two functions to define the extent and resolution of costmap grid if needed
function get_grid_extent()
    return 10;
end
function get_grid_resolution()
    return 0.1
end











