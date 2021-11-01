#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <tf/tf.h>
#include <Eigen/Dense>

geometry_msgs::Twist twist_sub_msg;

void twist_sub_cb(const geometry_msgs::Twist::ConstPtr& twist)
{
	twist_sub_msg = *twist;
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
	ros::init(argc, argv, "twist_conv");
	ros::NodeHandle nh("~");

	twist_sub_msg = geometry_msgs::Twist();

	ros::Subscriber sub = nh.subscribe<geometry_msgs::Twist>(nh.param<std::string>("subscriber/topic", "/twist"), nh.param("subscriber/queue", 1000), twist_sub_cb);
	ros::Publisher pub = nh.advertise<geometry_msgs::Twist>(nh.param<std::string>("publisher/topic", "/cmd_vel"), nh.param("publisher/queue", 1000));

	Eigen::MatrixXd angular_transform_matrix = matrixparam(nh, "angular_transform_matrix", Eigen::MatrixXd(4, 4));
	Eigen::MatrixXd linear_transform_matrix = matrixparam(nh, "linear_transform_matrix", Eigen::MatrixXd(4, 4));

	ros::Rate rate(nh.param("rate", 50));

	while (ros::ok()) {
		geometry_msgs::Twist twist_pub_msg;

		rate.sleep();

		ros::spinOnce();

		twist_pub_msg = twist_sub_msg;

		Eigen::Vector4d angular_sub(twist_sub_msg.angular.x, twist_sub_msg.angular.y, twist_sub_msg.angular.z, 1.0);
		Eigen::Vector4d angular_pub = angular_transform_matrix*angular_sub;
		twist_pub_msg.angular.x = angular_pub(0);
		twist_pub_msg.angular.y = angular_pub(1);
		twist_pub_msg.angular.z = angular_pub(2);

		Eigen::Vector4d linear_sub(twist_sub_msg.linear.x, twist_sub_msg.linear.y, twist_sub_msg.linear.z, 1.0);
		Eigen::Vector4d linear_pub = linear_transform_matrix*linear_sub;
		twist_pub_msg.linear.x = linear_pub(0);
		twist_pub_msg.linear.y = linear_pub(1);
		twist_pub_msg.linear.z = linear_pub(2);

		pub.publish(twist_pub_msg);
	}

	return EXIT_SUCCESS;
}
