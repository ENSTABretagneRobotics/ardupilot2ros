#include <rclcpp/rclcpp.hpp>
#include <geometry_msgs/msg/twist_stamped.hpp>
#include <geometry_msgs/msg/vector3_stamped.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <tf2/LinearMath/Quaternion.hpp>
#include <tf2/LinearMath/Matrix3x3.hpp>
#include <tf2/LinearMath/Vector3.hpp>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>
#include <cmath>

int main(int argc, char** argv)
{
	rclcpp::init(argc, argv);
	auto nh = std::make_shared<rclcpp::Node>("velocity_conv");
	auto pub = nh->create_publisher<geometry_msgs::msg::TwistStamped>(
		nh->declare_parameter<std::string>("publisher.topic", "/vel"), 10);
	const auto topic = nh->declare_parameter<std::string>("subscriber.topic", "/navsat/vel");
	const auto type = nh->declare_parameter<std::string>("subscriber.type", "vector3");
	rclcpp::Subscription<geometry_msgs::msg::Vector3Stamped>::SharedPtr vector_sub;
	rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub;
	if (type == "vector3")
	{
		vector_sub = nh->create_subscription<geometry_msgs::msg::Vector3Stamped>(topic, rclcpp::SensorDataQoS(),
			[pub](geometry_msgs::msg::Vector3Stamped::ConstSharedPtr msg) {
				geometry_msgs::msg::TwistStamped out;
				out.header = msg->header;
				out.twist.linear = msg->vector;
				pub->publish(out);
			});
	}
	else if (type == "odometry")
	{
		odom_sub = nh->create_subscription<nav_msgs::msg::Odometry>(topic, rclcpp::SensorDataQoS(),
			[nh, pub](nav_msgs::msg::Odometry::ConstSharedPtr msg) {
				tf2::Quaternion q;
				tf2::fromMsg(msg->pose.pose.orientation, q);
				if (!std::isfinite(q.length2()) || q.length2() < 1e-12)
				{
					RCLCPP_WARN_THROTTLE(nh->get_logger(), *nh->get_clock(), 5000, "Odometry orientation is invalid");
					return;
				}
				q.normalize();
				tf2::Vector3 v(msg->twist.twist.linear.x, msg->twist.twist.linear.y, msg->twist.twist.linear.z);
				v = tf2::Matrix3x3(q) * v;
				geometry_msgs::msg::TwistStamped out;
				out.header = msg->header;
				out.twist.linear.x = v.x();
				out.twist.linear.y = v.y();
				out.twist.linear.z = v.z();
				pub->publish(out);
			});
	}
	else
		throw std::invalid_argument("subscriber.type must be vector3 or odometry");
	rclcpp::spin(nh);
	rclcpp::shutdown();
	return 0;
}
