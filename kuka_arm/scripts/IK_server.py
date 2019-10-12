#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

##########################################
###############  هام جداٌ #################
## الرجاء حذف جميع التعليقات العربية في حال ##
####### الرغبة في إستخدام هذا الكود ########
#########################################

####################################
#######  إستيراد الوحدات #######

import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import

#########################################
######### تعريف رموز معايير دي إتش #######
# d1:8 أطوال الأضلاغ
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')

# a0:7 إنزياح الأضلاغ
a0, a1, a2, a3, a4, a5, a6, = symbols('a0:7')

# alpha0:7 زاوية الإلتواء
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

#q1:q8 زوايا المفاصل
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')

##########################################
########## بناء جدول معايير دي إتش #########

DH_TABLE = {alpha0:     0,  a0:      0,  d1:  0.75,  q1:       q1,
			alpha1: -pi/2,  a1:   0.35,  d2:     0,  q2: -pi/2+q2,
			alpha2:     0,  a2:   1.25,  d3:     0,  q3:       q3,
			alpha3: -pi/2,  a3: -0.054,  d4:   1.5,  q4:       q4,
			alpha4:  pi/2,  a4:      0,  d5:     0,  q5:       q5,
			alpha5: -pi/2,  a5:      0,  d6:     0,  q6:       q6,
			alpha6:     0,  a6:      0,  d7: 0.303,  q7:        0}

####################################################
#########  قانون جيب التمام لمثلث معلوم الأضلاغ #########

# بالمعادلة (a,b,c) يعطى القانون لمثلث أضلاعه
#  c^2 = a^2 + b^2 - 2a*b*cos(theta)
cos_law = lambda a, b, c: (a**2 + b**2 - c**2) / (2 * a * b)

####################################################
######### خساب قوس جيب تمام الزاوية #######

#إستخدام متطابقات الدوال المثلثية العكسية
# لحساب فوس جيب تمام الزاوية من قوس ظل تمام الزاوية
cos_inv = lambda x: atan2(sqrt(1 - x**2), x)

####################################################
########## حساب مصفوفات الدوران ##################

def URDF2DH(r, p ,y):
    URDF2DH_ROT_X = Matrix([
        [      1,      0,       0],
        [      0, cos(r), -sin(r)],
        [      0, sin(r),  cos(r)]])

    URDF2DH_ROT_y = Matrix([
        [ cos(p),      0,  sin(p)],
        [      0,      1,       0],
        [-sin(p),      0,  cos(p)]])

    URDF2DH_ROT_z = Matrix([
        [ cos(r),-sin(r),       0],
        [ sin(r), cos(r),       0],
        [      0,      0,       1]])

    return URDF2DH_ROT_X , URDF2DH_ROT_y , URDF2DH_ROT_z

#######################################################
##########   مصفوفة التحويل المتجانس ##########

def TF_MATRIX(alpha, a, theta, q):
    TF_MAT = Matrix([
        [           cos(q),           -sin(q),           0,             a],
        [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
        [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
        [                0,                 0,           0,             1]])

    return TF_MAT

######################################################
########### حساب معادلات الحركة الأمامية  ################

def forward_kinematics():

    ###################################################
    #########  تعويض قيم معايير دي إتش  ################
    #########    لحساب مصفوفات الحركة  ################

    T0_1 = TF_MATRIX(alpha0, a0, d1, q1).subs(DH_TABLE)
	T1_2 = TF_MATRIX(alpha1, a1, d2, q2).subs(DH_TABLE)
	T2_3 = TF_MATRIX(alpha2, a2, d3, q3).subs(DH_TABLE)
	T3_4 = TF_MATRIX(alpha3, a3, d4, q4).subs(DH_TABLE)
	T4_5 = TF_MATRIX(alpha4, a4, d5, q5).subs(DH_TABLE)
	T5_6 = TF_MATRIX(alpha5, a5, d6, q6).subs(DH_TABLE)
	T6_grip = TF_MATRIX(alpha6, a6, d7, q7).subs(DH_TABLE)

    ###################################################
    ############ التحويلات المتجانسة المُركبة ##############

	T0_2 = T0_1 * T1_2 				# القاعدة للضلع 2
	T0_3 = T0_2 * T2_3 				# القاعدة للضلع 3
	T0_4 = T0_3 * T3_4 				# القاعدة للضلع 4
	T0_5 = T0_4 * T4_5 				# القاعدة للضلع 5
	T0_6 = T0_5 * T5_6 				# القاعدة للضلع 6

    ##################################################
    ################ مصفوفات الدوران ##################

    R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
    R3_6 = T3_4[0:3, 0:3] * T4_5[0:3, 0:3] * T5_6[0:3, 0:3]

    #################################################
    ##### تحويل مقبض اليد بالنسبة للقاعدة ###############

    T0_grip = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_grip

    # R3_0 ألحركة العكسية

    R3_0 = R0_3.transpose()

    return R3_0, T0_grip


#######################################################
############# معادلات الحركة العكسية #####################

def inverse_kinematics(pos, ori, R3_0):

    ##################################################
    ############### DH إلى URDF التحويل من ###########

    # الحصول على موقع و دوران المقبض من المحاكاة
    grip = Matrix(pos)
    ori  = Matrix(ori)
    ###############################################
    ##### تحويل إحداثيات مقبض اليد إلى دي إتش ########

    r, p ,y = symbols('r p y')
    rot_x, rot_y, rot_z = URDF2DH(r, p, y)

    ################################################
    #### تصحيح دوران مقبض اليد نسبة إحداثيات دي إتش ###

    URDF2DH_grip_rot_correction = rot_z.subs(y, pi) * rot_y.subs(p, -pi/2)

    #########################################
    ##### تصحيح الدوران عندما يكون مقبض اليد ####
    #####      في أي وضعية عشوائية       ####

    rot_grip = rot_Z * rot_y * rot_x
    rot_grip = rot_grip * URDF2DH_grip_rot_correction
    rot_grip = rot_grip.subs({'r': ori[0], 'p': ori[1], 'y': ori[2]})

    R0_6 = rot_grip

    grip2WC_Translation = Matrix([
                                 [0],
                                 [0],
                                 [DH_TABLE[d7]]])

    WC = grip - R0_6*grip2WC_Translation

    ##################################################
    ################ حساب المسافات و الزوايا ############

    # حساب إسقاط مركز المعصم على المستوى  الإحداثي للقاعده
    #ناقص إنزياح الضلع من المفصل الثاني x0 هي العنصر في إتحاه xc
    xc = sqrt(WC[0]**2 + WC[1]**2 - DH_TABLE[a1])
    # ناقص طول الضلع من المفصل الثاني y0 هي العنصر في إتحاه yc
    yc = WC[2] - DH_TABLE[d1]

    # حساب المسافة بين المفاصل بإعتبار المفصلين 4 و 6 متصلين بالمفصل 5
    d2_3 = DH_TABLE[a2]
    d3_5 = sqrt(DH_TABLE[d4]**2 + DH_TABLE[a3]**2)
    d2_5 = sqrt(xc**2 + yc**2)

    alpha = atan2(yc, xc)
    beta = abs(atan2(DH_TABLE[a3], DH_TABLE[d4]))

    cos_a = cos_law(d2_5, d2_3, d3_5)
    cos_b = cos_law(d2_3, d3_5, d2_5)
    cos_c = cos_law(d2_5, d3_5, d2_3)

    angle_a = cos_inv(cos_a)
    angle_b = cos_inv(cos_b)
    angle_c = cos_inv(cos_c)

    ####################################################
    ################ حساب ثيتا1 و ثيتا2 و ثيتا3 ############

    theta1 = atan2(WC[1],WC[2]).evalf()
    theta2 = (pi/2 - (angle_a + alpha)).evalf()
    theta3 = (pi/2 - (angle_b + beta )).evalf()

    ###################################################
    ########### R3_6  حساب مصفوقة الدوران ##############
    R3_0 = R3_0.evalf(subs={q1: theta1, q2:theta2, q3:theta3})
    R3_6 = R3_0 * rot_grip # مصفوفة الدوران من مفصل 3 إلى المقبض

    ###################################################
    ######### حساب ثيتا4 و ثيتا5 و ثينا6 ##################
	theta5 = atan2( sqrt(R3_6[0,2]**2 + R3_6[2,2]**2), R3_6[1,2] ).evalf()
	if (theta5 > pi) :
		theta4 = atan2(-R3_6[2,2], R3_6[0,2]).evalf()
		theta6 = atan2(R3_6[1,1], -R3_6[1,0]).evalf()
	else:
		theta4 = atan2(R3_6[2,2], -R3_6[0,2]).evalf()
		theta6 = atan2(-R3_6[1,1], R3_6[1,0]).evalf()

	print("theta1: ", theta1)
	print("theta2: ", theta2)
	print("theta3: ", theta3)
	print("theta4: ", theta4)
	print("theta5: ", theta5)
	print("theta6: ", theta6)

	return [theta1, theta2, theta3, theta4, theta5, theta6]

###########################
###########################

def handle_calculate_IK(req):
	rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
	if len(req.poses) < 1:
		print "No valid poses received"
		return -1
	else:
		RMSE_EE_MAT = np.zeros((len(req.poses),3))
		#إجراء عملية تحليل معادلات الحركة الأمامية
		R3_0, T0_EE = forward_kinematics()

        # تهيئة الرد من السيرفر
		joint_trajectory_list = []

		for i in xrange(0, len(req.poses)):
			joint_trajectory_point = JointTrajectoryPoint()

            ###################################################
            ####### إستخراج موقع و دوران مقبض اليد من المحاكاة ####

            # موقع مقبض اليد
            px = req.poses[i].position.x
			py = req.poses[i].position.y
			pz = req.poses[i].position.z
			# زوايا الدورن رباعية الأبعاد
			oa = req.poses[i].orientation.x
			ob = req.poses[i].orientation.y
			oc = req.poses[i].orientation.z
			od = req.poses[i].orientation.w
			# تحويل الزوايا من رباعية الإبعاد إلى زوايا يولر
			#  (roll, pitch, yaw) = دوران مقبض اليد
			(roll, pitch, yaw) = tf.transformations.euler_from_quaternion([oa, ob, oc, od])
			pos = [px, py, pz]
			ori = [roll, pitch, yaw]

            # IK request تعبئة الرد من
			joint_trajectory_point = JointTrajectoryPoint()
            # حساب الثيتا
			[theta1, theta2, theta3, theta4, theta5, theta6] = inverse_kinematics(pos, ori, R3_0)

			joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
			#joint_trajectory_point.positions = inverse_kinematics(pos, ori, R3_0)
			joint_trajectory_list.append(joint_trajectory_point)

		rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
		print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time))

		return CalculateIKResponse(joint_trajectory_list)


def IK_server():
	global start_time
	start_time = time()
	print ("\nStart time is %04.4f seconds" % start_time)

	# initialize node and declare calculate_ik service
	rospy.init_node('IK_server')
	s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
	print "Ready to receive an IK request"
	rospy.spin()

if __name__ == "__main__":
	IK_server()
