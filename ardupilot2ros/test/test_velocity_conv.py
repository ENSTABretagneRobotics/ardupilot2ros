import math
import os
import signal
import subprocess
import time

from ament_index_python.packages import get_package_prefix
from geometry_msgs.msg import TwistStamped, Vector3Stamped
from nav_msgs.msg import Odometry
import pytest
import rclpy


@pytest.mark.parametrize('source', ['vector3', 'odometry'])
def test_velocity_conversion(source):
    executable = os.path.join(
        get_package_prefix('ardupilot2ros'), 'lib', 'ardupilot2ros', 'velocity_conv_node')
    process = subprocess.Popen([
        executable, '--ros-args', '-p', 'subscriber.type:=' + source,
        '-p', 'subscriber.topic:=/velocity_test_input',
        '-p', 'publisher.topic:=/velocity_test_output',
    ])
    rclpy.init()
    node = rclpy.create_node('velocity_conversion_test')
    received = []
    subscription = node.create_subscription(
        TwistStamped, '/velocity_test_output', received.append, 10)
    msg_type = Vector3Stamped if source == 'vector3' else Odometry
    publisher = node.create_publisher(msg_type, '/velocity_test_input', 10)
    msg = msg_type()
    if source == 'vector3':
        msg.header.frame_id = 'world'
        msg.vector.x = 2.0
        msg.vector.y = -3.0
    else:
        msg.header.frame_id = 'odom'
        msg.pose.pose.orientation.z = math.sin(math.pi / 4)
        msg.pose.pose.orientation.w = math.cos(math.pi / 4)
        msg.twist.twist.linear.x = 2.0
    try:
        if source == 'odometry':
            deadline = time.monotonic() + 10
            while publisher.get_subscription_count() == 0 and time.monotonic() < deadline:
                rclpy.spin_once(node, timeout_sec=0.1)
            assert publisher.get_subscription_count() > 0
            invalid = Odometry()
            invalid.pose.pose.orientation.w = 0.0
            invalid.twist.twist.linear.x = 2.0
            for _ in range(5):
                publisher.publish(invalid)
                rclpy.spin_once(node, timeout_sec=0.1)
            assert not received
        deadline = time.monotonic() + 15
        while not received and time.monotonic() < deadline:
            assert process.poll() is None
            publisher.publish(msg)
            rclpy.spin_once(node, timeout_sec=0.1)
        assert received
        result = received[-1]
        assert result.header.frame_id == msg.header.frame_id
        if source == 'vector3':
            assert result.twist.linear.x == pytest.approx(2.0)
            assert result.twist.linear.y == pytest.approx(-3.0)
        else:
            assert result.twist.linear.x == pytest.approx(0.0, abs=1e-6)
            assert result.twist.linear.y == pytest.approx(2.0)
        assert result.twist.angular.z == 0.0
    finally:
        node.destroy_subscription(subscription)
        node.destroy_node()
        rclpy.shutdown()
        process.send_signal(signal.SIGINT)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
