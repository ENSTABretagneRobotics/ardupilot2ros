#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <tf/tf.h>
#include <Eigen/Dense>

sensor_msgs::Imu imu_sub_msg;

void imu_sub_cb(const sensor_msgs::Imu::ConstPtr& imu)
{
	imu_sub_msg = *imu;
}

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

int main(int argc, char** argv) {
	ros::init(argc, argv, "imu_conv");
	ros::NodeHandle nh("~");

	imu_sub_msg = sensor_msgs::Imu();

	ros::Subscriber sub = nh.subscribe<sensor_msgs::Imu>(nh.param<std::string>("subscriber/topic", "/imu"), nh.param("subscriber/queue", 1000), imu_sub_cb);
	ros::Publisher pub = nh.advertise<sensor_msgs::Imu>(nh.param<std::string>("publisher/topic", "/imu/data"), nh.param("publisher/queue", 1000));

	Eigen::MatrixXd angular_velocity_transform_matrix = matrixparam(nh, "angular_velocity_transform_matrix", Eigen::MatrixXd(4, 4));
	Eigen::MatrixXd linear_acceleration_transform_matrix = matrixparam(nh, "linear_acceleration_transform_matrix", Eigen::MatrixXd(4, 4));
	Eigen::MatrixXd RPY_transform_matrix = matrixparam(nh, "RPY_transform_matrix", Eigen::MatrixXd(4, 4));

	ros::Rate rate(nh.param("rate", 50));

	while (ros::ok()) {
		sensor_msgs::Imu imu_pub_msg;

		rate.sleep();

		ros::spinOnce();

		imu_pub_msg = imu_sub_msg;

		Eigen::Vector4d angular_velocity_sub(imu_sub_msg.angular_velocity.x, imu_sub_msg.angular_velocity.y, imu_sub_msg.angular_velocity.z, 1.0);
		Eigen::Vector4d angular_velocity_pub = angular_velocity_transform_matrix*angular_velocity_sub;
		imu_pub_msg.angular_velocity.x = angular_velocity_pub(0);
		imu_pub_msg.angular_velocity.y = angular_velocity_pub(1);
		imu_pub_msg.angular_velocity.z = angular_velocity_pub(2);

		Eigen::Vector4d linear_acceleration_sub(imu_sub_msg.linear_acceleration.x, imu_sub_msg.linear_acceleration.y, imu_sub_msg.linear_acceleration.z, 1.0);
		Eigen::Vector4d linear_acceleration_pub = linear_acceleration_transform_matrix*linear_acceleration_sub;
		imu_pub_msg.linear_acceleration.x = linear_acceleration_pub(0);
		imu_pub_msg.linear_acceleration.y = linear_acceleration_pub(1);
		imu_pub_msg.linear_acceleration.z = linear_acceleration_pub(2);

		tf::Quaternion q;
		double roll = 0, pitch = 0, yaw = 0;
		tf::quaternionMsgToTF(imu_sub_msg.orientation, q);
		tf::Matrix3x3(q).getRPY(roll, pitch, yaw);
		Eigen::Vector4d RPY_sub(roll, pitch, yaw, 1.0);
		Eigen::Vector4d RPY_pub = RPY_transform_matrix*RPY_sub;
		roll = RPY_pub(0); pitch = RPY_pub(1); yaw = RPY_pub(2);
		q.setRPY(roll, pitch, yaw);
		imu_pub_msg.orientation.x = q[0];
		imu_pub_msg.orientation.y = q[1];
		imu_pub_msg.orientation.z = q[2];
		imu_pub_msg.orientation.w = q[3];

		pub.publish(imu_pub_msg);
	}

	return EXIT_SUCCESS;
}
