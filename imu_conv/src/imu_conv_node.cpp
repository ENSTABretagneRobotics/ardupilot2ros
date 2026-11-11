#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <tf2/LinearMath/Quaternion.hpp>
#include <tf2/LinearMath/Matrix3x3.hpp>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>
#include <Eigen/Dense>

sensor_msgs::msg::Imu imu_sub_msg;

void imu_sub_cb(const sensor_msgs::msg::Imu::ConstSharedPtr& imu)
{
	imu_sub_msg = *imu;
}

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

int main(int argc, char** argv) {
	rclcpp::init(argc, argv);
	auto nh = std::make_shared<rclcpp::Node>("imu_conv");

	imu_sub_msg = sensor_msgs::msg::Imu();
	imu_sub_msg.orientation.w = 1.0;

	auto sub = nh->create_subscription<sensor_msgs::msg::Imu>(nh->declare_parameter<std::string>("subscriber.topic", "/imu"), rclcpp::SensorDataQoS().keep_last(nh->declare_parameter("subscriber.queue", 1000)), imu_sub_cb);
	auto pub = nh->create_publisher<sensor_msgs::msg::Imu>(nh->declare_parameter<std::string>("publisher.topic", "/imu/data"), nh->declare_parameter("publisher.queue", 1000));

	Eigen::MatrixXd angular_velocity_transform_matrix = matrixparam(nh, "angular_velocity_transform_matrix", Eigen::MatrixXd::Identity(4, 4));
	Eigen::MatrixXd linear_acceleration_transform_matrix = matrixparam(nh, "linear_acceleration_transform_matrix", Eigen::MatrixXd::Identity(4, 4));
	Eigen::MatrixXd RPY_transform_matrix = matrixparam(nh, "RPY_transform_matrix", Eigen::MatrixXd::Identity(4, 4));

	rclcpp::WallRate rate(nh->declare_parameter("rate", 50));
	rclcpp::executors::SingleThreadedExecutor executor;
	executor.add_node(nh);

	while (rclcpp::ok()) {
		sensor_msgs::msg::Imu imu_pub_msg;

		rate.sleep();

		if (!rclcpp::ok()) break;
		executor.spin_some();

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

		tf2::Quaternion q;
		double roll = 0, pitch = 0, yaw = 0;
		tf2::fromMsg(imu_sub_msg.orientation, q);
		tf2::Matrix3x3(q).getRPY(roll, pitch, yaw);
		Eigen::Vector4d RPY_sub(roll, pitch, yaw, 1.0);
		Eigen::Vector4d RPY_pub = RPY_transform_matrix*RPY_sub;
		roll = RPY_pub(0); pitch = RPY_pub(1); yaw = RPY_pub(2);
		q.setRPY(roll, pitch, yaw);
		imu_pub_msg.orientation.x = q[0];
		imu_pub_msg.orientation.y = q[1];
		imu_pub_msg.orientation.z = q[2];
		imu_pub_msg.orientation.w = q[3];

		pub->publish(imu_pub_msg);
	}

	rclcpp::shutdown();
	return EXIT_SUCCESS;
}
