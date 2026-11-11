from launch import LaunchDescription
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory
import os


def generate_launch_description():
    config = os.path.join(get_package_share_directory('ardupilot2ros'), 'config', 'saturne.yaml')
    return LaunchDescription([
        Node(package='ardupilot2ros', executable='ardupilot2ros_node',
             name='ardupilot2ros', parameters=[config], output='screen', respawn=True),
    ])
