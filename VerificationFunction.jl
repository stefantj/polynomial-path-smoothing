
#Function verifyActuateablePath checks if the smooth path is feasible given robots limits
#Assumptions
# Robot follows the quadrotor model set out by Minimum  Snap  Trajectory  Generation  and  Control  for  Quadrotors - Mellinger
# Euler angles of Z-X-Y
# Polynomials of the same degree
#Inputs
# solution - an object containing points, and times
# max_vel - the maximum velocity that the robot is limited too in ros
# max_motor_rpm - the maximum rpm that a motor can get
#Outputs
# did_pass - a boolean that is true when the path passes
# timeProbv - a vector of times where motion is infeasible based on velocity
# timeProbm - a vector of times where motion is infeasible based on velocity
function verifyActuateablePath(solution::PolySol, max_vel::Float64, max_motor_rpm::Float64)
    #Extract important information from the solution object
    degree = 2 + 2*solution.params.cont_order #create the degree of each polynomial assuming 2 pts for each
    num_poly = solution.num_segs;
    xcoeffs = solution.x_coeffs;
    ycoeffs = solution.y_coeffs;
    zcoeffs = solution.z_coeffs;
    pcoeffs = solution.p_coeffs;
    time_vec = solution.times;
    timeProbv = [];
    timeProbm = [];
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
    timeRep = xdd;
    #Needed initializations
    z_B = 0; #variable to hold the up body vector of copter
    #Needed for coefficients and angular acceleration
    yawdd = xdd;
    yaw = xdd;

    #Loop through all segments
    for looper = 1:(num_poly)
        #Create time vector
        t = linspace(time_vec[looper],time_vec[looper+1],time_res);
        #update coeffs for ever poly
        xddcoeffs = xcoeffs[(red_degree_by+1:degree)+degree*(looper-1)];
        yddcoeffs = ycoeffs[(red_degree_by+1:degree)+degree*(looper-1)];
        zddcoeffs = zcoeffs[(red_degree_by+1:degree)+degree*(looper-1)];
        pddcoeffs = pcoeffs[(red_degree_by+1:degree)+degree*(looper-1)];
        
        for loop=1:time_res
            cleaner = (looper-1)*time_res+loop;
            timeRep[cleaner] = t[loop];
            #Calculate the velocity in positions
            for p=(2:degree)+degree*(looper-1)
                xd[cleaner] += xcoeffs[p]*t[loop]^(p-2-degree*(looper-1))*(p-1-degree*(looper-1));
                yd[cleaner] += ycoeffs[p]*t[loop]^(p-2-degree*(looper-1))*(p-1-degree*(looper-1));
                zd[cleaner] += zcoeffs[p]*t[loop]^(p-2-degree*(looper-1))*(p-1-degree*(looper-1));
            end
            #do a check if the max speed has been exceeded
            if( max_vel < sqrt((xd[cleaner])^2 + (yd[cleaner])^2 + (zd[cleaner])^2))
                did_pass = false;
                timeProbv = [timeProb t[loop]];
            end
            #Calculate the accelerations at every point
            for p=(3:degree)+degree*(looper-1)
                xdd[cleaner] += xcoeffs[p]*t[loop]^(p-3-degree*(looper-1))*(p-1-degree*(looper-1))*(p-2-degree*(looper-1));
                ydd[cleaner] += ycoeffs[p]*t[loop]^(p-3-degree*(looper-1))*(p-1-degree*(looper-1))*(p-2-degree*(looper-1));
                zdd[cleaner] += zcoeffs[p]*t[loop]^(p-3-degree*(looper-1))*(p-1-degree*(looper-1))*(p-2-degree*(looper-1));
                yawdd[cleaner] += pcoeffs[p]*t[loop]^(p-3-degree*(looper-1))*(p-1-degree*(looper-1))*(p-2-degree*(looper-1));
            end
            #Calculate the jerks in position at every point
            for p=(2:red_degree)
                xddd[cleaner] += xddcoeffs[p]*t[loop]^(p-2)*(p-1)*(p)*(1+p);
                yddd[cleaner] += yddcoeffs[p]*t[loop]^(p-2)*(p-1)*(p)*(1+p);
                zddd[cleaner] += zddcoeffs[p]*t[loop]^(p-2)*(p-1)*(p)*(1+p);
            end
            #Calculate the yaw at every point
            for p=(1:degree)+degree*(looper-1)
                yaw[cleaner] += pcoeffs[p]*t[loop]^(p-(1+degree*(looper-1)));
            end
            #Calculate the body centered vector pointing up relative to robot body
            z_B = [xdd[cleaner],ydd[cleaner],zdd[cleaner]]/
                sqrt((xdd[cleaner])^2 + (ydd[cleaner])^2 + (zdd[cleaner]+gravity)^2);
            #Calculate the x_c vector
            x_c = [sin(yaw[cleaner]),cos(yaw[cleaner]),0.0];
            #calculate y_B
            y_B = cross(z_B,x_c);
            #calculate x_B
            x_B = cross(y_B, z_B)
            #Calculate u1
            u1 = mass*sqrt((xdd[cleaner])^2 + (ydd[cleaner])^2 + (zdd[cleaner]+gravity)^2);
            #Calculate a_dot
            a_dot = [xddd[cleaner],yddd[cleaner],zddd[cleaner]];
            #Calculate h_w
            h_w = mass/u1*(a_dot-dot(z_B,a_dot)*z_B);
            #Calculate w_bc
            w_bc = -dot(h_w,y_B)*x_B + dot(h_w,x_B)*y_B + yaw[cleaner]*dot(z_w,z_B)*z_B;
            u2u3u4 = I*yawdd[cleaner]*z_w+cross(w_bc,I*w_bc);
            u1_vec[cleaner] = u1;
            u2_vec[cleaner] = u2u3u4[1];
            u3_vec[cleaner] = u2u3u4[2];
            u4_vec[cleaner] = u2u3u4[3];      

        end

    end

    #Calculate the needed rpms
    u_vec =[u1_vec'
        u2_vec'
        u3_vec'
        u4_vec'];
    rpms = u2rpms*u_vec
    println(size(rpms,1))
    println(size(rpms,2))
    figure();
    plot((1:num_poly*time_res),rpms[1,:])
    #check that rpms are not above a certain threshold
    for looper = 1:4
        timeProbm = [timeProbm timeRep[find(rpms[p,:].>max_motor_rpm)]]
    end
    
    #If timeProb v and m are not empty say so, velocity problems are to be notified first then motor probs
    if(!isempty(timeProbv))
        println("The path will break the speed limit")
        return timeProbv;
    end
    if(!isempty(timeProbm))
        println("The path will break the motor limits") 
        return timeProbv;
    end
    return [];
    
end


