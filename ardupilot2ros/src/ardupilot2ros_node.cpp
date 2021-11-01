#include "MAVLinkDevice.h"
#include "VectorNavInterface.h"
#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/NavSatFix.h>
#include <tf/tf.h>
#include <Eigen/Dense>

sensor_msgs::Imu imu_msg;
sensor_msgs::NavSatFix fix_msg;
geometry_msgs::TwistStamped vel_msg;

// Inspired from https://github.com/cra-ros-pkg/robot_localization/blob/32896d6d1aaec5a92c3a65e5f16dfc0b859a7d26/src/ros_filter.cpp#L1556 and 
// https://github.com/cra-ros-pkg/robot_localization/blob/32896d6d1aaec5a92c3a65e5f16dfc0b859a7d26/params/ukf_template.yaml#L209
// Should check more parsing step...?
Eigen::MatrixXd matrixparam(ros::NodeHandle nh, const std::string param, Eigen::MatrixXd matrixdefault)
{
	Eigen::MatrixXd matrix(matrixdefault);
	XmlRpc::XmlRpcValue matrix_tmp;

	nh.getParam(param, matrix_tmp);

	if (matrix_tmp.getType() != XmlRpc::XmlRpcValue::TypeArray) return matrixdefault;

	for (int i = 0; i < matrix.rows(); i++)
	{
		for (int j = 0; j < matrix.cols(); j++)
		{
			// These matrices can cause problems if all the types
			// aren't specified with decimal points. Handle that
			// using string streams.
			std::ostringstream ostr;
			ostr << matrix_tmp[matrix.cols()*i+j];
			std::istringstream istr(ostr.str());
			istr >> matrix(i, j);
		}
	}

	return matrix;
}

void imuCb(const sensor_msgs::Imu::ConstPtr& imu)
{
	//ROS_INFO("I heard: [%f]", imu->angular_velocity.z);
	imu_msg = *imu;
}

void fixCb(const sensor_msgs::NavSatFix::ConstPtr& fix)
{
	//ROS_INFO("I heard: [%f] [%f]", fix->latitude, fix->longitude);
	fix_msg = *fix;
}

void velCb(const geometry_msgs::TwistStamped::ConstPtr& vel)
{
	//ROS_INFO("I heard: [%f] [%f] [%f]", vel->twist.linear.x, vel->twist.linear.y, vel->twist.linear.z);
	vel_msg = *vel;
}

int main(int argc, char** argv) {

	THREAD_IDENTIFIER MAVLinkDeviceThreadId;
	THREAD_IDENTIFIER VectorNavInterfaceThreadId;

	ros::init(argc, argv, "ardupilot2ros");
	ros::NodeHandle nh("~");

	imu_msg = sensor_msgs::Imu();
	fix_msg = sensor_msgs::NavSatFix();
	vel_msg = geometry_msgs::TwistStamped();

	Eigen::MatrixXd cmd_vel_servos_to_angular_linear_transform_matrix = matrixparam(nh, "cmd_vel/servos_to_angular_linear_transform_matrix", Eigen::MatrixXd(7, 9));

	Eigen::MatrixXd imu_angular_velocity_transform_matrix = matrixparam(nh, "imu/angular_velocity_transform_matrix", Eigen::MatrixXd(4, 4));
	Eigen::MatrixXd imu_linear_acceleration_transform_matrix = matrixparam(nh, "imu/linear_acceleration_transform_matrix", Eigen::MatrixXd(4, 4));
	Eigen::MatrixXd imu_RPY_transform_matrix = matrixparam(nh, "imu/RPY_transform_matrix", Eigen::MatrixXd(4, 4));

	Eigen::MatrixXd vel_angular_transform_matrix = matrixparam(nh, "vel/angular_transform_matrix", Eigen::MatrixXd(4, 4));
	Eigen::MatrixXd vel_linear_transform_matrix = matrixparam(nh, "vel/linear_transform_matrix", Eigen::MatrixXd(4, 4));

	bool bOverrideFixStatus(nh.param("bOverrideFixStatus", true));

	// VectorNavInterface and MAVLinkDevice files are currently shared with another project, so there are some internal parameters to set and unnessary code...

	strcpy(szVectorNavInterfacePath, (nh.param<std::string>("szVectorNavInterfacePath", "127.0.0.1:5765")).c_str());
	VectorNavInterfaceBaudRate = (nh.param("VectorNavInterfaceBaudRate", 230400));
	VectorNavInterfaceTimeout = (nh.param("VectorNavInterfaceTimeout", 2000));

	// Might not be necessary to change those parameters depending on the area since the conversions made will be in both directions,
	// so errors related to the fact Earth is not flat should be cancelled...?
	lat_env = (nh.param("lat_env", 48.418079));
	long_env = (nh.param("long_env", -4.473487));
	alt_env = (nh.param("alt_env", 88.0));
	angle_env = (M_PI/2.0-nh.param("angle_env", 90.0)*M_PI/180.0);
	AirPressure = (nh.param("AirPressure", 1.0));

	robid = (nh.param("robid", BUGGY_ROBID)); // Normally not used if bExternal is 1 in MAVLinkDevice0.txt...

	InitGlobals();

	ros::Publisher cmd_vel_pub = nh.advertise<geometry_msgs::Twist>(nh.param<std::string>("cmd_vel/topic", "/cmd_vel"), nh.param("cmd_vel/queue", 1000));
	ros::Subscriber imu_sub = nh.subscribe<sensor_msgs::Imu>(nh.param<std::string>("imu/topic", "/imu/data"), nh.param("imu/queue", 1000), imuCb);
	ros::Subscriber fix_sub = nh.subscribe<sensor_msgs::NavSatFix>(nh.param<std::string>("fix/topic", "/fix"), nh.param("fix/queue", 1000), fixCb);
	ros::Subscriber vel_sub = nh.subscribe<geometry_msgs::TwistStamped>(nh.param<std::string>("vel/topic", "/vel"), nh.param("vel/queue", 1000), velCb);

	ros::Duration(2.0).sleep();

	CreateDefaultThread(MAVLinkDeviceThread, NULL, &MAVLinkDeviceThreadId);
	DetachThread(MAVLinkDeviceThreadId); // Not easy to stop it correctly...
	CreateDefaultThread(VectorNavInterfaceThread, NULL, &VectorNavInterfaceThreadId);
	DetachThread(VectorNavInterfaceThreadId); // Not easy to stop it correctly...

	ros::Rate rate(nh.param("rate", 50));

	while (ros::ok()) {
		sensor_msgs::Imu imu_enu_msg;
		geometry_msgs::TwistStamped vel_enu_msg;
		geometry_msgs::Twist cmd_vel_msg;

		rate.sleep();

		ros::spinOnce();

		EnterCriticalSection(&StateVariablesCS);

		// imu_msg, fix_sub_msg, vel_sub_msg to VectorNav..

		// Conversions to ENU coordinate system...
		Eigen::Vector4d vel_angular = vel_angular_transform_matrix*Eigen::Vector4d(vel_msg.twist.angular.x, vel_msg.twist.angular.y, vel_msg.twist.angular.z, 1.0);
		vel_enu_msg.twist.angular.x = vel_angular(0); vel_enu_msg.twist.angular.y = vel_angular(1); vel_enu_msg.twist.angular.z = vel_angular(2);
		Eigen::Vector4d vel_linear = vel_linear_transform_matrix*Eigen::Vector4d(vel_msg.twist.linear.x, vel_msg.twist.linear.y, vel_msg.twist.linear.z, 1.0);
		vel_enu_msg.twist.linear.x = vel_linear(0); vel_enu_msg.twist.linear.y = vel_linear(1); vel_enu_msg.twist.linear.z = vel_linear(2);
		Eigen::Vector4d imu_angular_velocity = imu_angular_velocity_transform_matrix*Eigen::Vector4d(imu_msg.angular_velocity.x, imu_msg.angular_velocity.y, imu_msg.angular_velocity.z, 1.0);
		imu_enu_msg.angular_velocity.x = imu_angular_velocity(0); imu_enu_msg.angular_velocity.y = imu_angular_velocity(1); imu_enu_msg.angular_velocity.z = imu_angular_velocity(2);
		Eigen::Vector4d imu_linear_acceleration = imu_linear_acceleration_transform_matrix*Eigen::Vector4d(imu_msg.linear_acceleration.x, imu_msg.linear_acceleration.y, imu_msg.linear_acceleration.z, 1.0);
		imu_enu_msg.linear_acceleration.x = imu_linear_acceleration(0); imu_enu_msg.linear_acceleration.y = imu_linear_acceleration(1); imu_enu_msg.linear_acceleration.z = imu_linear_acceleration(2);
		tf::Quaternion q;
		double roll = 0, pitch = 0, yaw = 0;
		tf::quaternionMsgToTF(imu_msg.orientation, q);
		tf::Matrix3x3(q).getRPY(roll, pitch, yaw);
		Eigen::Vector4d imu_RPY = imu_RPY_transform_matrix*Eigen::Vector4d(roll, pitch, yaw, 1.0);
		roll = imu_RPY(0); pitch = imu_RPY(1); yaw = imu_RPY(2);
		q.setRPY(roll, pitch, yaw);
		imu_enu_msg.orientation.x = q[0]; imu_enu_msg.orientation.y = q[1]; imu_enu_msg.orientation.z = q[2]; imu_enu_msg.orientation.w = q[3];

		// Intermediate variables used by VectorNavInterface...
		double x = 0, y = 0, z = 0;
		GPS2EnvCoordSystem(lat_env, long_env, alt_env, angle_env, fix_msg.latitude, fix_msg.longitude, fix_msg.altitude, &x, &y, &z);
		x_gps = x; y_gps = y; z_gps = z;
		cog_gps = fmod_2PI(0*M_PI/2.0+atan2(vel_enu_msg.twist.linear.y, vel_enu_msg.twist.linear.x)-angle_env);
		sog = sqrt(pow(vel_enu_msg.twist.linear.x, 2)+pow(vel_enu_msg.twist.linear.y, 2));
		vx_ned = vel_enu_msg.twist.linear.x; vy_ned = vel_enu_msg.twist.linear.y; vz_ned = vel_enu_msg.twist.linear.z;
		xhat = x; yhat = y; zhat = z;
		if (bOverrideFixStatus)
		{
			GNSSqualitySimulator = AUTONOMOUS_GNSS_FIX;
			GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 20;
			GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 1.0;
		}
		else
		{
			switch (fix_msg.status.status)
			{
			case sensor_msgs::NavSatStatus::STATUS_FIX:
			case sensor_msgs::NavSatStatus::STATUS_SBAS_FIX:
				GNSSqualitySimulator = AUTONOMOUS_GNSS_FIX;
				GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 20;
				GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 1.0;
				break;
			case sensor_msgs::NavSatStatus::STATUS_GBAS_FIX:
				GNSSqualitySimulator = RTK_FLOAT;
				GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 20;
				GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 0.8;
				break;
			case sensor_msgs::NavSatStatus::STATUS_NO_FIX:
			default:
				GNSSqualitySimulator = GNSS_NO_FIX;
				GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 0;
				GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 0;
				break;
			}
		}
		//tf::Quaternion q;
		//double roll = 0, pitch = 0, yaw = 0;
		//tf::quaternionMsgToTF(imu_enu_msg.orientation, q);
		//tf::Matrix3x3(q).getRPY(roll, pitch, yaw);
		phihat = fmod_2PI(roll); thetahat = fmod_2PI(pitch); psihat = fmod_2PI(0*M_PI/2.0+yaw-angle_env);
		omegaxhat = imu_enu_msg.angular_velocity.x; omegayhat = imu_enu_msg.angular_velocity.y; omegazhat = imu_enu_msg.angular_velocity.z;
		accrxhat = imu_enu_msg.linear_acceleration.x; accryhat = imu_enu_msg.linear_acceleration.y; accrzhat = imu_enu_msg.linear_acceleration.z;
		//ROS_INFO("[%f] [%f] [%f] [%f] [%f] [%f] [%f]", x, y, z, roll, pitch, yaw, cog_gps);

		// Output...
		Eigen::MatrixXd cmd_vel_servos = Eigen::MatrixXd(9, 1);
		cmd_vel_servos(0) = u1_servo_out_MAVLinkDevice[0]; cmd_vel_servos(1) = u2_servo_out_MAVLinkDevice[0]; cmd_vel_servos(2) = u3_servo_out_MAVLinkDevice[0];
		cmd_vel_servos(3) = u4_servo_out_MAVLinkDevice[0]; cmd_vel_servos(4) = u5_servo_out_MAVLinkDevice[0]; cmd_vel_servos(5) = u6_servo_out_MAVLinkDevice[0];
		cmd_vel_servos(6) = u7_servo_out_MAVLinkDevice[0]; cmd_vel_servos(7) = u8_servo_out_MAVLinkDevice[0]; cmd_vel_servos(8) = 1.0;
		Eigen::MatrixXd cmd_vel_angular_linear = cmd_vel_servos_to_angular_linear_transform_matrix*cmd_vel_servos;
		cmd_vel_msg.angular.x = cmd_vel_angular_linear(0); cmd_vel_msg.angular.y = cmd_vel_angular_linear(1); cmd_vel_msg.angular.z = cmd_vel_angular_linear(2);
		cmd_vel_msg.linear.x = cmd_vel_angular_linear(3); cmd_vel_msg.linear.y = cmd_vel_angular_linear(4); cmd_vel_msg.linear.z = cmd_vel_angular_linear(5);
		//ROS_INFO("[%f] [%f] [%f] [%f] [%f] [%f] [%f] [%f] [%f]", 
		//	cmd_vel_servos(0), cmd_vel_servos(1), cmd_vel_servos(2), cmd_vel_servos(3), cmd_vel_servos(4), cmd_vel_servos(5), cmd_vel_servos(6), cmd_vel_servos(7), cmd_vel_servos(8));
		//ROS_INFO("[%f] [%f] [%f] [%f] [%f] [%f] [%f]", 
		//	cmd_vel_angular_linear(0), cmd_vel_angular_linear(1), cmd_vel_angular_linear(2), cmd_vel_angular_linear(3), cmd_vel_angular_linear(4), cmd_vel_angular_linear(5), cmd_vel_angular_linear(6));
		//std::cout << cmd_vel_servos << endl;
		//std::cout << endl;
		//std::cout << cmd_vel_servos_to_angular_linear_transform_matrix << endl;
		//std::cout << endl;
		//std::cout << cmd_vel_angular_linear << endl;
		//std::cout << endl;

		LeaveCriticalSection(&StateVariablesCS);

		cmd_vel_pub.publish(cmd_vel_msg);
	}

	ReleaseGlobals();

	return EXIT_SUCCESS;
}
