import os

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def generate_launch_description():
    namespace = LaunchConfiguration('robot_namespace')
    prefix = ['/', namespace]
    config = os.path.join(get_package_share_directory('ardupilot2ros'), 'config', 'warthog_cpr_sim.yaml')
    return LaunchDescription([
        DeclareLaunchArgument('robot_namespace', default_value='w200_0000'),
        Node(package='ardupilot2ros', executable='ardupilot2ros_node',
             name='ardupilot2ros', output='screen', respawn=True,
             parameters=[config, {
                 'use_sim_time': True,
                 'cmd_vel.topic': [*prefix, '/cmd_vel'],
                 'imu.topic': [*prefix, '/sensors/imu_0/data'],
                 'fix.topic': [*prefix, '/sensors/gps_0/fix'],
             }]),
        Node(package='ardupilot2ros', executable='velocity_conv_node',
             name='velocity_conv', output='screen', respawn=True,
             parameters=[{
                 'use_sim_time': True,
                 'subscriber.type': 'odometry',
                 'subscriber.topic': [*prefix, '/platform/odom'],
                 'publisher.topic': '/vel',
             }]),
    ])
