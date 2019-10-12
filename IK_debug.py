from sympy import *
from time import time
from mpmath import radians, degrees
import tf

'''
Format of test case is [ [ [EE position], [EE orientation as quaternions] ], [WC location], [joint angles] ]
You can generate additional test cases by setting up your kuka project and running 
`$ roslaunch kuka_arm forward_kinematics.launch`
From here you can adjust the joint angles to find thetas, use the gripper to extract 
positions and orientation (in quaternion xyzw) and lastly use link 5 to find the position of the 
wrist center. These newly generated test cases can be added to the test_cases dictionary.
'''

test_cases = {1:[[[2.16135,-1.42635,1.55109],
                  [0.708611,0.186356,-0.157931,0.661967]],
                  [1.89451,-1.44302,1.69366],
                  [-0.65,0.45,-0.36,0.95,0.79,0.49]],
              2:[[[-0.56754,0.93663,3.0038],
                  [0.62073, 0.48318,0.38759,0.480629]],
                  [-0.638,0.64198,2.9988],
                  [-0.79,-0.11,-2.33,1.94,1.14,-3.68]],
              3:[[[-1.3863,0.02074,0.90986],
                  [0.01735,-0.2179,0.9025,0.371016]],
                  [-1.1669,-0.17989,0.85137],
                  [-2.99,-0.12,0.94,4.06,1.29,-4.12]],
              4:[],
              5:[]}


def test_code(test_case):
    ## Set up code
    ## Do not modify!
    x = 0
    class Position:
        def __init__(self,EE_pos):
            self.x = EE_pos[0]
            self.y = EE_pos[1]
            self.z = EE_pos[2]
    class Orientation:
        def __init__(self,EE_ori):
            self.x = EE_ori[0]
            self.y = EE_ori[1]
            self.z = EE_ori[2]
            self.w = EE_ori[3]

    position = Position(test_case[0][0])
    orientation = Orientation(test_case[0][1])

    class Combine:
        def __init__(self,position,orientation):
            self.position = position
            self.orientation = orientation

    comb = Combine(position,orientation)

    class Pose:
        def __init__(self,comb):
            self.poses = [comb]

    req = Pose(comb)
    start_time = time()
    
    ########################################################################################
    ## 

    ########################################################################################
    ########################           FORWARD KINEMATICS           ########################
    ################          >> GIVEN JOINT ANGLES OF ALL JOINTS           ################
    ##########          >> FIND END-EFFECTOR POSITION AND ORIENTATION             ##########
    ########################################################################################
    
    ########################################################
    ############## DEFINE DH PARAMETER SYMBOLS #############
    ########################################################
    # Link lengtha: d1:8
    d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
    # Link offseta: a0:7
    a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
    # Twist angles: alpha0:7
    alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
    # Joint angles: q1:8
    q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')

    ########################################################
    ######### CREATE DH PARAMETER TABLE/DICTIONARY #########
    ########################################################
    # DH Parameter for KUKA210 Robot Manipulator
    DH_TABLE = {alpha0:     0,  a0:      0,  d1:  0.75,  q1:       q1, 
                alpha1: -pi/2,  a1:   0.35,  d2:     0,  q2: -pi/2+q2,
                alpha2:     0,  a2:   1.25,  d3:     0,  q3:       q3,
                alpha3: -pi/2,  a3: -0.054,  d4:   1.5,  q4:       q4,
                alpha4:  pi/2,  a4:      0,  d5:     0,  q5:       q5,
                alpha5: -pi/2,  a5:      0,  d6:     0,  q6:       q6,
                alpha6:     0,  a6:      0,  d7: 0.303,  q7:        0}

    ########################################################
    ####### DEFINE HOMOGENEOUS TRANSFORMATION MATRIX #######
    ########################################################
    def TF_MATRIX(alpha, a, d, q):
        TF_MAT = Matrix([[            cos(q),               -sin(q),               0,                a],
                         [ sin(q)*cos(alpha),     cos(q)*cos(alpha),     -sin(alpha),    -sin(alpha)*d],
                         [ sin(q)*sin(alpha),     cos(q)*sin(alpha),      cos(alpha),     cos(alpha)*d],
                         [                 0,                     0,               0,                1]])
        return TF_MAT

    '''
    # Base Link to Link-1
    T0_1 = Matrix([[            cos(q1),               -sin(q1),               0,                 0],
                   [sin(q1)*cos(alpha0),    cos(q1)*cos(alpha0),    -sin(alpha0),   -sin(alpha0)*d1],
                   [sin(q1)*sin(alpha0),    cos(q1)*sin(alpha0),     cos(alpha0),    cos(alpha0)*d1],
                   [                  0,                      0,               0,                 1]])
    # Link-1 to Link-2
    T1_2 = Matrix([[            cos(q2),               -sin(q2),               0,                 0],
                   [sin(q2)*cos(alpha1),    cos(q2)*cos(alpha1),    -sin(alpha1),   -sin(alpha1)*d2],
                   [sin(q2)*sin(alpha1),    cos(q2)*sin(alpha1),     cos(alpha1),    cos(alpha1)*d2],
                   [                  0,                      0,               0,                 1]])
    # Link-2 to Link-3
    T2_3 = Matrix([[            cos(q3),               -sin(q3),               0,                 0],
                   [sin(q3)*cos(alpha2),    cos(q3)*cos(alpha2),    -sin(alpha2),   -sin(alpha2)*d3],
                   [sin(q3)*sin(alpha2),    cos(q3)*sin(alpha2),     cos(alpha2),    cos(alpha2)*d3],
                   [                  0,                      0,               0,                 1]])
    # Link-3 to Link-4
    T3_4 = Matrix([[            cos(q4),               -sin(q4),               0,                 0],
                   [sin(q4)*cos(alpha3),    cos(q4)*cos(alpha3),    -sin(alpha3),   -sin(alpha3)*d4],
                   [sin(q4)*sin(alpha3),    cos(q4)*sin(alpha3),     cos(alpha3),    cos(alpha3)*d4]
                   [                  0,                      0,               0,                 1]])
    # Link-4 to Link-5
    T4_5 = Matrix([[            cos(q5),               -sin(q5),               0,                 0],
                   [sin(q5)*cos(alpha4),    cos(q5)*cos(alpha4),    -sin(alpha4),   -sin(alpha4)*d5],
                   [sin(q5)*sin(alpha4),    cos(q5)*sin(alpha4),     cos(alpha4),    cos(alpha4)*d5],
                   [                  0,                      0,               0,                 1]])
    # Link-5 to Link-6
    T5_6 = Matrix([[            cos(q6),               -sin(q6),               0,                 0],
                   [sin(q6)*cos(alpha5),    cos(q6)*cos(alpha5),    -sin(alpha5),   -sin(alpha5)*d6],
                   [sin(q6)*sin(alpha5),    cos(q6)*sin(alpha5),     cos(alpha5),    cos(alpha5)*d6],
                   [                  0,                      0,               0,                 1]])
    # Link-6 to End Effector
    T6_EE = Matrix([[            cos(q7),               -sin(q7),               0,                 0],
                    [sin(q7)*cos(alpha6),    cos(q7)*cos(alpha6),    -sin(alpha6),   -sin(alpha6)*d7],
                    [sin(q7)*sin(alpha6),    cos(q7)*sin(alpha6),     cos(alpha6),    cos(alpha6)*d7],
                    [                  0,                      0,               0,                 1]])
    '''
    ########################################################
    ###### EVALUATE ALL TRANSFORMATION MATRIXCES AND #######
    ###### SUBSTITUTE PARAMETER VALUES FROM DH-TABLE #######
    ########################################################
    T0_1 = TF_MATRIX(alpha0, a0, d1, q1).subs(DH_TABLE)
    T1_2 = TF_MATRIX(alpha1, a1, d2, q2).subs(DH_TABLE)
    T2_3 = TF_MATRIX(alpha2, a2, d3, q3).subs(DH_TABLE)
    T3_4 = TF_MATRIX(alpha3, a3, d4, q4).subs(DH_TABLE)
    T4_5 = TF_MATRIX(alpha4, a4, d5, q5).subs(DH_TABLE)
    T5_6 = TF_MATRIX(alpha5, a5, d6, q6).subs(DH_TABLE)
    T6_EE = TF_MATRIX(alpha6, a6, d7, q7).subs(DH_TABLE)

    ########################################################
    ######## COMPOSITION OF HOMOGENEOUS TRANSFORMS #########
    ########################################################
    T0_2 = simplify(T0_1 * T1_2) # base_link to link-2
    T0_3 = simplify(T0_2 * T2_3) # base_link to link-3
    T0_4 = simplify(T0_3 * T3_4) # base_link to link-4
    T0_5 = simplify(T0_4 * T4_5) # base_link to link-5
    T0_6 = simplify(T0_5 * T5_6) # base_link to link-6
    
    ########################################################
    ################## ROTATION MATRICES ###################
    ########################################################
    R0_1 = T0_1[0:3, 0:3]        # base_link to link-1
    R0_2 = R0_1 * T1_2[0:3, 0:3] # base_link to link-2
    R0_3 = R0_2 * T2_3[0:3, 0:3] # base_link to link-3
    R0_4 = R0_3 * T3_4[0:3, 0:3] # base_link to link-4
    R0_5 = R0_4 * T4_5[0:3, 0:3] # base_link to link-5
    R0_6 = R0_5 * T5_6[0:3, 0:3] # base_link to link-6

    RR3_6 = T3_4[0:3, 0:3]*T4_5[0:3, 0:3]*T5_6[0:3, 0:3]  #R3_4*R4_5*R5_6
    print("R3_4 = ", T3_4[0:3, 0:3])
    print("R4_5 = ", T4_5[0:3, 0:3])
    print("R5_6 = ", T5_6[0:3, 0:3])
    print("R3_6 = ", RR3_6)
    
    ########################################################
    ####### END-EFFECTOR TRANSFORM W.R.T. BASE-LINK ########
    ########################################################
    T0_EE = simplify(T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_EE)
    '''
    # Print matrices
    print("T0_1 = ", T0_1.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("T0_2 = ", T0_2.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("T0_3 = ", T0_3.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("T0_4 = ", T0_4.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("T0_5 = ", T0_5.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("T0_6 = ", T0_6.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("T0_EE = ", T0_EE.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))

    print("R0_1 = ", R0_1.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("R0_2 = ", R0_2.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("R0_3 = ", R0_3.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("R0_4 = ", R0_4.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("R0_5 = ", R0_5.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    print("R0_6 = ", R0_1.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
    '''
    
    ########################################################################################
    ####################               INVERSE KINEMATICS             ######################
    ##########          >> GIVEN END-EFFECTOR POSITION AND ORIENTATION      ################
    ##########                >> FIND JOINT ANGLES OF ALL JOINTS            ################
    ########################################################################################
    ## Insert IK code here!

    theta1 = 0
    theta2 = 0
    theta3 = 0
    theta4 = 0
    theta5 = 0
    theta6 = 0

    '''
    Since the last three joints in our robot are revolute and their joint axes intersect 
    at a single point, we have a case of spherical wrist with joint_5 being the common 
    intersection point and hence the wrist center (WC). This allows us to kinematically 
    decouple the IK problem into Inverse Position and Inverse Orientation problems.

    First let us solve for the Inverse Position problem. Since we have the case of a 
    spherical wrist involving joints 4,5,6, the position of the wrist center is governed 
    by the first three joints. We can obtain the position of the wrist center by using the 
    complete transformation matrix we derived based on the end-effector pose.
    '''

    ########################################################
    ########## USE ATAN2 TO GET SIGNED THETA ###############
    ########## FROM X=COS(THETA)  ##########################
    ########################################################
    cos_inv = lambda X: atan2(sqrt(1 - X**2), X)
    #cos_inv = lambda X: acos(X)

    
    ########################################################
    ################## LAW OF COSINES ######################
    ########################################################
    # Law of cosines for any given triangle: c^2 = a^2 + b^2 - 2a*b*cos(theta)
    cos_law = lambda a, b, c: (a**2 + b**2 - c**2) / (2 * a * b) # Returns cosine of angle between sides a and b


    ########################################################
    ########## EXTRACT END-EFFECTOR POSITION AND ###########
    ########## ORIENTATION FROM SIMULATION REQUEST #########
    ########################################################
    px = req.poses[x].position.x
    py = req.poses[x].position.y
    pz = req.poses[x].position.z
    # Orientation angles in quaternions
    oa = req.poses[x].orientation.x
    ob = req.poses[x].orientation.y
    oc = req.poses[x].orientation.z
    od = req.poses[x].orientation.w
    # px = position.x
    # py = position.y
    # pz = position.z
    # # Orientation angles in quaternions
    # oa = orientation.x
    # ob = orientation.y
    # oc = orientation.z
    # od = orientation.w
    # Convert orientation angles from quaternions to Euler angles
    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion([oa, ob, oc, od])

    ########################################################
    ###### GRIPPER/EE COORDNIATE TRANSFORMATION FROM #######
    ###### URDF COORDINATE FRAME TO DH COORDINATE FRAME ####
    ########################################################
    # Gripper/End-Effector Rotation Matrices [URDF -> DH]
    r, p, y = symbols('r p y')
    URDF2DH_ROT_X = Matrix([[      1,           0,          0 ],
                            [      0,      cos(r),    -sin(r) ],
                            [      0,      sin(r),     cos(r) ]])

    URDF2DH_ROT_Y = Matrix([[ cos(p),           0,     sin(p) ],
                            [      0,           1,          0 ],
                            [-sin(p),           0,     cos(p) ]])

    URDF2DH_ROT_Z = Matrix([[ cos(y),     -sin(y),          0 ],
                            [ sin(y),      cos(y),          0 ],
                            [      0,           0,          1 ]])

    ########################################################
    ###### BODY-FIXED ROTATION OF GRIPPER/EE [UDRF >> DH] ##
    ########################################################   
    URDF2DH_EE_ROT_CORRECTION = URDF2DH_ROT_Z.subs(y, pi) * URDF2DH_ROT_Y.subs(p, -pi/2)
    '''
    URDF2DH_ROT_Z = Matrix([[-1,      0,          0 ],
                            [ 0,     -1,          0 ],
                            [ 0,      0,          1 ]])

    URDF2DH_ROT_Y = Matrix([[ 0,     0,          -1 ],
                            [ 0,     1,           0 ],
                            [ 1,     0,           0 ]])
    URDF2DH_EE_ROT_CORRECTION = [[0, 0, 1]
                                 [0,-1, 0]
                                 [1, 0, 0]]
    '''
    ########################################################
    ###### FOR ARBITRARY EE ORIENTATION FROM SIMULATION ####
    ###### CONVERT GRIPPER/EE ORIENTATION [UDRF >> DH] #####
    ########################################################
    URDF2DH_ROT_EE = URDF2DH_ROT_Z * URDF2DH_ROT_Y * URDF2DH_ROT_X
    URDF2DH_ROT_EE = URDF2DH_ROT_EE * URDF2DH_EE_ROT_CORRECTION
    URDF2DH_ROT_EE = URDF2DH_ROT_EE.subs({'r': roll, 'p': pitch, 'y': yaw})

    #T0_EE_URDF2DH = simplify(T0_EE * URDF2DH_EE_ROT_CORRECTION) 
    print("(px, py, pz) = ", (px, py, pz) )
    print("d7 = ", DH_TABLE[d7])
    print("URDF2DH_ROT_EE = ", URDF2DH_ROT_EE[:,2])


    ########################################################
    ###################### URDF TO DH ######################
    ########################################################
    EE = Matrix([[px],
                 [py],
                 [pz]])

    EE2WC_TRANSLATION = Matrix([[0],
                                [0],
                                [DH_TABLE[d7]]])

    R0_6 = URDF2DH_ROT_EE   # The Gripper link is rigidly connected to Link-6 >> q7=0 >> R6_EE = I
                            # URDF2DH_ROT_EE = R0_EE = R0_6*R6_EE = R0_6*I = R0_6

    #WC = EE - DH_TABLE[d7]*URDF2DH_ROT_EE[:,2] # Wrist Center w.r.t. Base Link
    WC = EE - R0_6*EE2WC_TRANSLATION # Wrist Center w.r.t. Base Link

    print("WC = ", WC)

    ########################################################
    ############# CALCULATE DISTANCE AND ANGLES ############
    ########################################################
    # Calculate wrist center on projected X0-Y0 plane (xc, yc) from base frame coordinates.
    # xc is the component in X0 direction minus offset from joint 2
    xc = sqrt(WC[0]**2 + WC[1]**2) - DH_TABLE[a1] #norm([w[0], w[1]]) - s[a1]
    # yc is the component in Z0 direction minus offset from joint 2
    yc = WC[2] - DH_TABLE[d1]

    # Calculate distances between joints
    d2_3 = DH_TABLE[a2]                               # distance between joint-2 to joint-3
    d3_5 = sqrt(DH_TABLE[a3]**2 + DH_TABLE[d4]**2)    # distance between joint-3 to joint-5/WC
    d2_5 = sqrt( xc**2 + yc**2 )                      # distance between joint-2 to joint-5/WC 
    
    alpha = atan2(yc, xc)
    beta = atan2(-DH_TABLE[a3], sqrt(DH_TABLE[d4]**2 - DH_TABLE[a3]**2))
    #beta = atan2(-DH_TABLE[a3], DH_TABLE[d4])
    
    '''
    cos_a = (d2_5*d2_5 + d2_3*d2_3 - d3_5*d3_5) / (2*d2_5*d2_3)
    cos_b = (d2_3*d2_3 + d3_5*d3_5 - d2_5*d2_5) / (2*d2_3*d3_5)
    cos_c = (d3_5*d3_5 + d2_5*d2_5 - d2_3*d2_3) / (2*d3_5*d2_5)
    '''
    cos_a = cos_law(d2_5, d2_3, d3_5)
    cos_b = cos_law(d2_3, d3_5, d2_5)
    cos_c = cos_law(d3_5, d2_5, d2_3)

    angle_a = cos_inv(cos_a)
    angle_b = cos_inv(cos_b)
    angle_c = cos_inv(cos_c)

    theta1 = atan2(WC[1], WC[0]).evalf()
    theta2 = (pi/2 - angle_a - alpha).evalf()
    theta3 = (pi/2 - angle_b - beta).evalf()

    R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

    R3_0 = R0_3.transpose()         # inverse(R0_3) = R3_0 = transpose(R0_3)

    R3_6 = R3_0 * URDF2DH_ROT_EE    # Rotation matrix from wrist to gripper
    print("R0_3 = ", R0_3)
    print("R3_6 = ", R3_6) 
    '''
    # R3_6 =   
    [[-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6), -sin(q4)*cos(q6) - cos(q4)*cos(q5)*sin(q6), -cos(q4)*sin(q5)],
     [                           sin(q5)*cos(q6),                           -sin(q5)*sin(q6),          cos(q5)],
     [-sin(q4)*cos(q5)*cos(q6) - cos(q4)*sin(q6),  sin(q4)*cos(q5)*sin(q6) - cos(q4)*cos(q6),  sin(q4)*sin(q5)]]

    R3_6 [theta5 = 0 or theta5 == 2*pi] =   
    [[ cos(q4)*cos(q6), -cos(q4)*sin(q6),   0],
     [               0,                0,   1],
     [-sin(q4)*cos(q6),  sin(q4)*sin(q6),   0]]

    R3_6 [theta5 = pi] =   
    [[-cos(q4)*cos(q6),  cos(q4)*sin(q6),   0],
     [               0,                0,  -1],
     [ sin(q4)*cos(q6), -sin(q4)*sin(q6),   0]]
    
    '''

    theta5 = atan2( sqrt(R3_6[0,2]**2 + R3_6[2,2]**2), R3_6[1,2] ).evalf() 
    if (theta5 > pi) :
        theta4 = atan2(-R3_6[2,2], R3_6[0,2]).evalf() 
        theta6 = atan2(R3_6[1,1],-R3_6[1,0]).evalf() 
    else:
        theta4 = atan2(R3_6[2,2], -R3_6[0,2]).evalf() 
        theta6 = atan2(-R3_6[1,1],R3_6[1,0]).evalf() 

    '''
    if (theta5 == 0 or theta5 == 2*pi):
        theta4 = atan2(-R3_6[2,0], R3_6[0,0]).evalf() 
        theta6 = atan2(-R3_6[0,1], R3_6[0,0]).evalf() 
    elif (theta5 == pi):
        theta4 = atan2(R3_6[2,0], -R3_6[0,0]).evalf() 
        theta6 = atan2(R3_6[0,1], -R3_6[0,0]).evalf() 
    elif ( pi < theta5 < 2*pi):
        theta4 = atan2(-R3_6[2,2], R3_6[0,2]).evalf() 
        theta6 = atan2(R3_6[1,1], -R3_6[1,0]).evalf() 
    else:
        theta4 = atan2(R3_6[2,2], -R3_6[0,2]).evalf() 
        theta6 = atan2(-R3_6[1,1], R3_6[1,0]).evalf() 
    '''
    print("theta4: ", theta4)
    print("theta5: ", theta5)
    print("theta6: ", theta6)

    FK = T0_EE.evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4: theta4, q5: theta5, q6: theta6})
    print("FK = ", FK)
    ## 
    ########################################################################################
    
    ########################################################################################
    ## For additional debugging add your forward kinematics here. Use your previously calculated thetas
    ## as the input and output the position of your end effector as your_ee = [x,y,z]

    ## (OPTIONAL) YOUR CODE HERE!


    ## End your code input for forward kinematics here!
    ########################################################################################

    ## For error analysis please set the following variables of your WC location and EE location in the format of [x,y,z]
    your_wc = [WC[0], WC[1], WC[2]] # <--- Load your calculated WC values in this array
    your_ee = [FK[0,3], FK[1,3], FK[2,3]] # <--- Load your calculated end effector value from your forward kinematics
    ########################################################################################

    ## Error analysis
    print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time))

    # Find WC error
    if not(sum(your_wc)==3):
        wc_x_e = abs(your_wc[0]-test_case[1][0])
        wc_y_e = abs(your_wc[1]-test_case[1][1])
        wc_z_e = abs(your_wc[2]-test_case[1][2])
        wc_offset = sqrt(wc_x_e**2 + wc_y_e**2 + wc_z_e**2)
        print ("\nWrist error for x position is: %04.8f" % wc_x_e)
        print ("Wrist error for y position is: %04.8f" % wc_y_e)
        print ("Wrist error for z position is: %04.8f" % wc_z_e)
        print ("Overall wrist offset is: %04.8f units" % wc_offset)

    # Find theta errors
    t_1_e = abs(theta1-test_case[2][0])
    t_2_e = abs(theta2-test_case[2][1])
    t_3_e = abs(theta3-test_case[2][2])
    t_4_e = abs(theta4-test_case[2][3])
    t_5_e = abs(theta5-test_case[2][4])
    t_6_e = abs(theta6-test_case[2][5])
    print ("\nTheta 1 error is: %04.8f" % t_1_e)
    print ("Theta 2 error is: %04.8f" % t_2_e)
    print ("Theta 3 error is: %04.8f" % t_3_e)
    print ("Theta 4 error is: %04.8f" % t_4_e)
    print ("Theta 5 error is: %04.8f" % t_5_e)
    print ("Theta 6 error is: %04.8f" % t_6_e)
    print ("\n**These theta errors may not be a correct representation of your code, due to the fact \
           \nthat the arm can have muliple positions. It is best to add your forward kinmeatics to \
           \nconfirm whether your code is working or not**")
    print (" ")

    # Find FK EE error
    if not(sum(your_ee)==3):
        ee_x_e = abs(your_ee[0]-test_case[0][0][0])
        ee_y_e = abs(your_ee[1]-test_case[0][0][1])
        ee_z_e = abs(your_ee[2]-test_case[0][0][2])
        ee_offset = sqrt(ee_x_e**2 + ee_y_e**2 + ee_z_e**2)
        print ("\nEnd effector error for x position is: %04.8f" % ee_x_e)
        print ("End effector error for y position is: %04.8f" % ee_y_e)
        print ("End effector error for z position is: %04.8f" % ee_z_e)
        print ("Overall end effector offset is: %04.8f units \n" % ee_offset)




if __name__ == "__main__":
    # Change test case number for different scenarios
    test_case_number = 3

    test_code(test_cases[test_case_number])
