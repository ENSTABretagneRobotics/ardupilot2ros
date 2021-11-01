#include "MAVLinkDevice.h"
#include "VectorNavInterface.h"
#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/NavSatFix.h>
#include "tf/tf.h"
#include "Eigen/Dense"

sensor_msgs::Imu imu_sub_msg;
sensor_msgs::NavSatFix fix_sub_msg;
geometry_msgs::TwistStamped vel_sub_msg;

void imuCb(const sensor_msgs::Imu::ConstPtr& imu)
{
    //ROS_INFO("I heard: [%f]", imu->angular_velocity.z);
	imu_sub_msg = *imu;
}

void fixCb(const sensor_msgs::NavSatFix::ConstPtr& fix)
{
    //ROS_INFO("I heard: [%f] [%f]", fix->latitude, fix->longitude);
	fix_sub_msg = *fix;
}

void velCb(const geometry_msgs::TwistStamped::ConstPtr& vel)
{
    //ROS_INFO("I heard: [%f] [%f] [%f]", vel->twist.linear.x, vel->twist.linear.y, vel->twist.linear.z);
	vel_sub_msg = *vel;
}

int main(int argc, char** argv) {

    THREAD_IDENTIFIER MAVLinkDeviceThreadId;
    THREAD_IDENTIFIER VectorNavInterfaceThreadId;

    ros::init(argc, argv, "ardupilot2ros");
	ros::NodeHandle nh("~");

    ros::Publisher cmd_vel_pub = nh.advertise<geometry_msgs::Twist>(nh.param<std::string>("cmd_vel/topic", "/cmd_vel"), nh.param("cmd_vel/queue", 1000));
    ros::Subscriber imu_sub = nh.subscribe<sensor_msgs::Imu>(nh.param<std::string>("imu/topic", "/imu/data"), nh.param("imu/queue", 1000), imuCb);
    ros::Subscriber fix_sub = nh.subscribe<sensor_msgs::NavSatFix>(nh.param<std::string>("fix/topic", "/navsat/fix"), nh.param("fix/queue", 1000), fixCb);
    ros::Subscriber vel_sub = nh.subscribe<geometry_msgs::TwistStamped>(nh.param<std::string>("vel/topic", "/vel"), nh.param("vel/queue", 1000), velCb);

	double max_linear_x_speed(nh.param("max_linear_x_speed", 3.0));
	double max_angular_z_speed(nh.param("max_angular_z_speed", 3.0));

	strcpy(szVectorNavInterfacePath, (nh.param<std::string>("szVectorNavInterfacePath", "192.168.0.16:5765")).c_str());
	VectorNavInterfaceBaudRate = (nh.param("VectorNavInterfaceBaudRate", 230400));
	VectorNavInterfaceTimeout = (nh.param("VectorNavInterfaceTimeout", 2000));

    // Might not be necessary to change those parameters depending on the area since the conversions made will be in both directions...
   	lat_env = (nh.param("lat_env", 48.418079));
   	long_env = (nh.param("long_env", -4.473487));
   	alt_env = (nh.param("alt_env", 88.0));
   	angle_env = (M_PI/2.0-nh.param("angle_env", 90.0)*M_PI/180.0);
    AirPressure = (nh.param("AirPressure", 1.0));

    robid = BUGGY_ROBID; // For MAVLinkDevice...
    //target_followme = MAVLINKDEVICE0_TARGET; // For MAVLinkDevice...
    InitGlobals();

    ros::Duration(2.0).sleep();

    CreateDefaultThread(MAVLinkDeviceThread, NULL, &MAVLinkDeviceThreadId);
	DetachThread(MAVLinkDeviceThreadId); // Not easy to stop it correctly...
	CreateDefaultThread(VectorNavInterfaceThread, NULL, &VectorNavInterfaceThreadId);
	DetachThread(VectorNavInterfaceThreadId); // Not easy to stop it correctly...

    ros::Rate rate(nh.param("rate", 50));

    while (ros::ok()) {
        geometry_msgs::Twist cmd_vel_msg;

        rate.sleep();

        ros::spinOnce();

		EnterCriticalSection(&StateVariablesCS);
		// imu_sub_msg, fix_sub_msg, vel_sub_msg (assumed to be in NWU coordinate system) to VectorNav..
		double x = 0, y = 0, z = 0;
		GPS2EnvCoordSystem(lat_env, long_env, alt_env, angle_env, fix_sub_msg.latitude, fix_sub_msg.longitude, fix_sub_msg.altitude, &x, &y, &z);
        x_gps = x; y_gps = y; z_gps = z;
        cog_gps = fmod_2PI(M_PI/2.0+atan2(vel_sub_msg.twist.linear.y,vel_sub_msg.twist.linear.x)-angle_env);
        sog = sqrt(pow(vel_sub_msg.twist.linear.x,2)+pow(vel_sub_msg.twist.linear.y,2));
        GNSSqualitySimulator = AUTONOMOUS_GNSS_FIX;//RTK_FIXED;
        GPS_high_acc_nbsat = GPS_med_acc_nbsat = GPS_low_acc_nbsat = 20;
        GPS_high_acc_HDOP = GPS_med_acc_HDOP = GPS_low_acc_HDOP = 0.8;
        xhat = x; yhat = y; zhat = z;
        tf::Quaternion q;
		double roll = 0, pitch = 0, yaw = 0;
		tf::quaternionMsgToTF(imu_sub_msg.orientation, q);
		tf::Matrix3x3(q).getRPY(roll, pitch, yaw);
        phihat = fmod_2PI(roll); thetahat = fmod_2PI(pitch); psihat = fmod_2PI(M_PI/2.0+yaw-angle_env);
        omegaxhat = imu_sub_msg.angular_velocity.x; omegayhat = imu_sub_msg.angular_velocity.y; omegazhat = imu_sub_msg.angular_velocity.z;
        accrxhat = imu_sub_msg.linear_acceleration.x; accryhat = imu_sub_msg.linear_acceleration.y; accrzhat = imu_sub_msg.linear_acceleration.z;
        //ROS_INFO("I heard: [%f] [%f] [%f] [%f] [%f] [%f]", x, y, z, roll, pitch, yaw);

        // Assume ArduRover...
		cmd_vel_msg.linear.x = max_linear_x_speed*u_servo_out_MAVLinkDevice[0];
		cmd_vel_msg.angular.z = max_angular_z_speed*uw_servo_out_MAVLinkDevice[0];
		LeaveCriticalSection(&StateVariablesCS);

        cmd_vel_pub.publish(cmd_vel_msg);
    }

	ReleaseGlobals();

	return EXIT_SUCCESS;
}
