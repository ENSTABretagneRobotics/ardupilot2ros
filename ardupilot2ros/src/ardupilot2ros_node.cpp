#include "MAVLinkDevice.h"
#include "VectorNavInterface.h"
#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <geometry_msgs/msg/twist_stamped.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/nav_sat_fix.hpp>
#include <tf2/LinearMath/Quaternion.hpp>
#include <tf2/LinearMath/Matrix3x3.hpp>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>
#include <Eigen/Dense>

sensor_msgs::msg::Imu imu_msg;
sensor_msgs::msg::NavSatFix fix_msg;
geometry_msgs::msg::TwistStamped vel_msg;

Eigen::MatrixXd matrixparam(const rclcpp::Node::SharedPtr& nh, const std::string& param, const Eigen::MatrixXd& matrixdefault)
{
	std::vector<double> defaults;
	for (int i = 0; i < matrixdefault.rows(); i++)
		for (int j = 0; j < matrixdefault.cols(); j++)
			defaults.push_back(matrixdefault(i, j));
	const auto values = nh->declare_parameter<std::vector<double>>(param, defaults);
	if (values.size() != defaults.size())
		throw std::invalid_argument(param + " must contain " + std::to_string(defaults.size()) + " values");
	Eigen::MatrixXd matrix(matrixdefault.rows(), matrixdefault.cols());
	for (int i = 0; i < matrix.rows(); i++)
		for (int j = 0; j < matrix.cols(); j++)
			matrix(i, j) = values[matrix.cols()*i+j];
	return matrix;
}

void imuCb(const sensor_msgs::msg::Imu::ConstSharedPtr& imu)
{
	//ROS_INFO("I heard: [%f]", imu->angular_velocity.z);
	imu_msg = *imu;
}

void fixCb(const sensor_msgs::msg::NavSatFix::ConstSharedPtr& fix)
{
	//ROS_INFO("I heard: [%f] [%f]", fix->latitude, fix->longitude);
	fix_msg = *fix;
}

void velCb(const geometry_msgs::msg::TwistStamped::ConstSharedPtr& vel)
{
	//ROS_INFO("I heard: [%f] [%f] [%f]", vel->twist.linear.x, vel->twist.linear.y, vel->twist.linear.z);
	vel_msg = *vel;
}

int main(int argc, char** argv) {

	THREAD_IDENTIFIER MAVLinkDeviceThreadId;
	THREAD_IDENTIFIER VectorNavInterfaceThreadId;

	rclcpp::init(argc, argv);
	auto nh = std::make_shared<rclcpp::Node>("ardupilot2ros");

	imu_msg = sensor_msgs::msg::Imu();
	imu_msg.orientation.w = 1.0;
	fix_msg = sensor_msgs::msg::NavSatFix();
	vel_msg = geometry_msgs::msg::TwistStamped();

	Eigen::MatrixXd cmd_vel_servos_to_angular_linear_transform_matrix = matrixparam(nh, "cmd_vel.servos_to_angular_linear_transform_matrix", Eigen::MatrixXd::Zero(7, 9));

	Eigen::MatrixXd imu_angular_velocity_transform_matrix = matrixparam(nh, "imu.angular_velocity_transform_matrix", Eigen::MatrixXd::Identity(4, 4));
	Eigen::MatrixXd imu_linear_acceleration_transform_matrix = matrixparam(nh, "imu.linear_acceleration_transform_matrix", Eigen::MatrixXd::Identity(4, 4));
	Eigen::MatrixXd imu_RPY_transform_matrix = matrixparam(nh, "imu.RPY_transform_matrix", Eigen::MatrixXd::Identity(4, 4));

	Eigen::MatrixXd vel_angular_transform_matrix = matrixparam(nh, "vel.angular_transform_matrix", Eigen::MatrixXd::Identity(4, 4));
	Eigen::MatrixXd vel_linear_transform_matrix = matrixparam(nh, "vel.linear_transform_matrix", Eigen::MatrixXd::Identity(4, 4));

	bool bOverrideFixStatus(nh->declare_parameter("bOverrideFixStatus", true));

	// VectorNavInterface and MAVLinkDevice files are currently shared with another project, so there are some internal parameters to set and unnessary code...

	strcpy(szVectorNavInterfacePath, (nh->declare_parameter<std::string>("szVectorNavInterfacePath", "127.0.0.1:5765")).c_str());
	VectorNavInterfaceBaudRate = (nh->declare_parameter("VectorNavInterfaceBaudRate", 230400));
	VectorNavInterfaceTimeout = (nh->declare_parameter("VectorNavInterfaceTimeout", 2000));

	// Might not be necessary to change those parameters depending on the area since the conversions made will be in both directions,
	// so errors related to the fact Earth is not flat should be cancelled...?
	lat_env = (nh->declare_parameter("lat_env", 48.418079));
	long_env = (nh->declare_parameter("long_env", -4.473487));
	alt_env = (nh->declare_parameter("alt_env", 88.0));
	angle_env = (M_PI/2.0-nh->declare_parameter("angle_env", 90.0)*M_PI/180.0);
	AirPressure = (nh->declare_parameter("AirPressure", 1.0));

	robid = (nh->declare_parameter("robid", BUGGY_ROBID)); // Normally not used if bExternal is 1 in MAVLinkDevice0.txt...

	InitGlobals();

	const auto cmd_vel_topic = nh->declare_parameter<std::string>("cmd_vel.topic", "/cmd_vel");
	const auto cmd_vel_queue = nh->declare_parameter("cmd_vel.queue", 1000);
	const bool cmd_vel_stamped = nh->declare_parameter("cmd_vel.stamped", false);
	rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel_pub;
	rclcpp::Publisher<geometry_msgs::msg::TwistStamped>::SharedPtr cmd_vel_stamped_pub;
	if (cmd_vel_stamped)
		cmd_vel_stamped_pub = nh->create_publisher<geometry_msgs::msg::TwistStamped>(cmd_vel_topic, cmd_vel_queue);
	else
		cmd_vel_pub = nh->create_publisher<geometry_msgs::msg::Twist>(cmd_vel_topic, cmd_vel_queue);
	auto imu_sub = nh->create_subscription<sensor_msgs::msg::Imu>(nh->declare_parameter<std::string>("imu.topic", "/imu/data"), rclcpp::SensorDataQoS().keep_last(nh->declare_parameter("imu.queue", 1000)), imuCb);
	auto fix_sub = nh->create_subscription<sensor_msgs::msg::NavSatFix>(nh->declare_parameter<std::string>("fix.topic", "/fix"), rclcpp::SensorDataQoS().keep_last(nh->declare_parameter("fix.queue", 1000)), fixCb);
	auto vel_sub = nh->create_subscription<geometry_msgs::msg::TwistStamped>(nh->declare_parameter<std::string>("vel.topic", "/vel"), rclcpp::SensorDataQoS().keep_last(nh->declare_parameter("vel.queue", 1000)), velCb);

	rclcpp::sleep_for(std::chrono::seconds(2));

	CreateDefaultThread(MAVLinkDeviceThread, NULL, &MAVLinkDeviceThreadId);
	DetachThread(MAVLinkDeviceThreadId); // Not easy to stop it correctly...
	CreateDefaultThread(VectorNavInterfaceThread, NULL, &VectorNavInterfaceThreadId);
	DetachThread(VectorNavInterfaceThreadId); // Not easy to stop it correctly...

	rclcpp::WallRate rate(nh->declare_parameter("rate", 50));
	rclcpp::executors::SingleThreadedExecutor executor;
	executor.add_node(nh);

	while (rclcpp::ok()) {
		sensor_msgs::msg::Imu imu_enu_msg;
		geometry_msgs::msg::TwistStamped vel_enu_msg;
		geometry_msgs::msg::Twist cmd_vel_msg;

		rate.sleep();

		if (!rclcpp::ok()) break;
		executor.spin_some();

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
		tf2::Quaternion q;
		double roll = 0, pitch = 0, yaw = 0;
		tf2::fromMsg(imu_msg.orientation, q);
		tf2::Matrix3x3(q).getRPY(roll, pitch, yaw);
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
			case sensor_msgs::msg::NavSatStatus::STATUS_FIX:
			case sensor_msgs::msg::NavSatStatus::STATUS_SBAS_FIX:
				GNSSqualitySimulator = AUTONOMOUS_GNSS_FIX;
				GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 20;
				GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 1.0;
				break;
			case sensor_msgs::msg::NavSatStatus::STATUS_GBAS_FIX:
				GNSSqualitySimulator = RTK_FLOAT;
				GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 20;
				GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 0.8;
				break;
			case sensor_msgs::msg::NavSatStatus::STATUS_NO_FIX:
			default:
				GNSSqualitySimulator = GNSS_NO_FIX;
				GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 0;
				GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 0;
				break;
			}
		}
		//tf2::Quaternion q;
		//double roll = 0, pitch = 0, yaw = 0;
		//tf2::fromMsg(imu_enu_msg.orientation, q);
		//tf2::Matrix3x3(q).getRPY(roll, pitch, yaw);
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

		if (cmd_vel_stamped)
		{
			geometry_msgs::msg::TwistStamped stamped;
			stamped.header.stamp = nh->now();
			stamped.twist = cmd_vel_msg;
			cmd_vel_stamped_pub->publish(stamped);
		}
		else
			cmd_vel_pub->publish(cmd_vel_msg);
	}

	ReleaseGlobals();

	rclcpp::shutdown();
	return EXIT_SUCCESS;
}
