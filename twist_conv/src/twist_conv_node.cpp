#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <geometry_msgs/msg/twist_stamped.hpp>
#include <Eigen/Dense>

geometry_msgs::msg::Twist twist_sub_msg;

void twist_sub_cb(const geometry_msgs::msg::Twist::ConstSharedPtr& twist)
{
	twist_sub_msg = *twist;
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
	auto nh = std::make_shared<rclcpp::Node>("twist_conv");

	twist_sub_msg = geometry_msgs::msg::Twist();

	auto sub = nh->create_subscription<geometry_msgs::msg::Twist>(nh->declare_parameter<std::string>("subscriber.topic", "/twist"), nh->declare_parameter("subscriber.queue", 1000), twist_sub_cb);
	auto pub = nh->create_publisher<geometry_msgs::msg::Twist>(nh->declare_parameter<std::string>("publisher.topic", "/cmd_vel"), nh->declare_parameter("publisher.queue", 1000));

	Eigen::MatrixXd angular_transform_matrix = matrixparam(nh, "angular_transform_matrix", Eigen::MatrixXd::Identity(4, 4));
	Eigen::MatrixXd linear_transform_matrix = matrixparam(nh, "linear_transform_matrix", Eigen::MatrixXd::Identity(4, 4));

	rclcpp::WallRate rate(nh->declare_parameter("rate", 50));
	rclcpp::executors::SingleThreadedExecutor executor;
	executor.add_node(nh);

	while (rclcpp::ok()) {
		geometry_msgs::msg::Twist twist_pub_msg;

		rate.sleep();

		if (!rclcpp::ok()) break;
		executor.spin_some();

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

		pub->publish(twist_pub_msg);
	}

	rclcpp::shutdown();
	return EXIT_SUCCESS;
}
