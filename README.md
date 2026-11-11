# ardupilot2ros

This workspace connects ROS sensor and velocity topics to an ArduPilot Rover
SITL instance through its VectorNav external AHRS and MAVLink interfaces. It
contains `ardupilot2ros`, `imu_conv`, and `twist_conv`.

## Install and build

Use Ubuntu 24.04 with the [binary ROS 2 Jazzy installation and ROS apt
repository](https://docs.ros.org/en/jazzy/Installation/Ubuntu-Install-Debians.html)
already set up. Then:

```bash
sudo apt update
sudo apt install -y git python3-colcon-common-extensions python3-rosdep \
  libeigen3-dev ros-jazzy-mavlink
# Run rosdep init once if this machine has not already been initialized.
test -e /etc/ros/rosdep/sources.list.d/20-default.list || sudo rosdep init
rosdep update
mkdir -p ~/ardupilot2ros_ws/src
cd ~/ardupilot2ros_ws/src
git clone -b dev https://github.com/ENSTABretagneRobotics/ardupilot2ros.git
cd ~/ardupilot2ros_ws
source /opt/ros/jazzy/setup.bash
rosdep install --from-paths src --ignore-src -r -y
colcon build --event-handlers console_direct+
source install/setup.bash
ros2 pkg list | grep -E '^(ardupilot2ros|imu_conv|twist_conv)$'
```

`ros-jazzy-mavlink` supplies the generated MAVLink v2 headers; a separate
MAVLink source checkout is not needed. `rosdep` installs the
remaining declared ROS dependencies, including tf2 and message packages.

## Bridge and configuration

The bridge consumes `sensor_msgs/msg/Imu`, `sensor_msgs/msg/NavSatFix`, and
`geometry_msgs/msg/TwistStamped` and publishes `geometry_msgs/msg/Twist` by
default. Set `cmd_vel.stamped` to true for a `TwistStamped` command instead.
Sensor subscriptions use sensor-data QoS. Parameters and transform matrices
are in `ardupilot2ros/config`; matrix values are row-major doubles. The
default matrices are identity for sensor transforms and zero for servo output,
so configure the servo matrix before driving a robot. The standalone
`imu_conv` and `twist_conv` nodes accept their installed YAML files through
`--ros-args --params-file`.

The MAVLink device reads `MAVLinkDevice0.txt` from its **current working
directory**, and writes logs under `log/` there. For example, after building:

```bash
cd ~/ardupilot2ros_ws
mkdir -p log
cp src/ardupilot2ros/ardupilot2ros/other/MAVLinkDevice0.txt .
cp src/ardupilot2ros/ardupilot2ros/other/eeprom.bin .
```

Inspect the connection settings in that file and the launch YAML. Run a
compatible ArduPilot Rover SITL in a separate terminal with the supplied
`eeprom.bin`, for example `./ardurover --home 48.418079,-4.473487,88,0
--model rover --speedup 1`. The historical parameters are in
`ardupilot2ros/other/ardupilot2ros.param`; the binary and SITL setup are not
provided by this repository. The default MAVLink and VectorNav connections
are `127.0.0.1:5760` and `127.0.0.1:5765`, respectively.

## Clearpath Warthog simulation

Install the [Clearpath Jazzy simulator](https://docs.clearpathrobotics.com/docs/ros/tutorials/simulator/install/)
and `ros-jazzy-ros-gz`, configure `~/clearpath/robot.yaml` for a W200 robot and
GPS, then start `ros2 launch clearpath_gz simulation.launch.py`. In another
terminal, from the workspace directory and after sourcing the workspace:

```bash
ros2 launch ardupilot2ros warthog_cpr_sim.launch.py robot_namespace:=w200_0000
```

The launch assumes `/<namespace>/cmd_vel` is `TwistStamped`,
`/<namespace>/sensors/imu_0/data` is `Imu`,
`/<namespace>/sensors/gps_0/fix` is `NavSatFix`, and
`/<namespace>/platform/odom` is `nav_msgs/msg/Odometry`. The GPS topic depends
on the configured sensor and generated Gazebo bridge. If the simulator
instead exposes `/<namespace>/sensors/gps_0/navsat`, bridge that Gazebo NavSat
topic with `ros_gz_bridge` and remap or change `fix.topic` to match.

The velocity converter rotates the odometry twist from its child frame to the
odometry frame using its pose orientation. This substitutes wheel odometry for
the old ROS 1 GPS velocity topic; it assumes the odometry world frame is ENU
and aligned with the simulated GPS world. Validate that assumption and the IMU
orientation, GPS status, QoS, command topic, and servo scaling for the actual
robot configuration before use. `velocity_conv_node` also supports the old
`Vector3Stamped` to `TwistStamped` conversion with `subscriber.type:=vector3`.

All three packages have been compiled in a clean Jazzy container. Standalone
conversion nodes were exercised with ROS messages, and the Warthog launch was
started without the simulator or SITL. The Clearpath Jazzy simulator generator
and control configuration were inspected for the interfaces above; the full
simulator, ArduPilot SITL and end-to-end Warthog operation have **not** been
runtime-tested. Simulator compatibility is therefore not established by
execution.
