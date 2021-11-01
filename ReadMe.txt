ardupilot2ros is a ROS node that enables the use of some ArduPilot functionality for ROS-based robots, without the need to change hardware.
It takes /imu/data, /fix and /vel typical ROS messages as inputs (robot angles, position, speed), transfer them to a simulated autopilot (see https://ardupilot.org/dev/docs/sitl-simulator-software-in-the-loop.html) using its VectorNav external AHRS interface (see https://github.com/ArduPilot/ardupilot/issues/11479, https://github.com/ArduPilot/ardupilot/blob/master/libraries/SITL/SIM_VectorNav.cpp), and then read the SERVO_OUTPUT_RAW MAVLink message from the simulated autopilot (PWM output values) and transfer it to the robot as /cmd_vel ROS message (robot inputs).

It can be tested with e.g. https://github.com/ENSTABretagneRobotics/SATURNE_simulation or https://clearpathrobotics.com/assets/guides/melodic/warthog/WarthogSimulation.html. http://wiki.ros.org/topic_tools/transform, imu_conv and twist_conv might be useful to make some conversions if needed.

Note that most of the code from this projet is currently shared with https://github.com/ENSTABretagneRobotics/UxVCtrl and https://github.com/ENSTABretagneRobotics/OSUtils so there is a lot of unused code. 
You might also need to put https://github.com/mavlink/c_library_v2 inside ardupilot2ros/src/mavlink.

Prepare 3 terminals. You might need to set ROS_MASTER_URI and ROS_IP accordingly if the terminals should not be on the same computer.

# First terminal
cd ~/catkin_ws
catkin_make
source devel/setup.bash
# For VMware virtual machines : export SVGA_VGPU10=0
roslaunch warthog_gazebo empty_world.launch
# For SATURNE_simulation : roslaunch mybot_gazebo mybot_world.launch

# Second terminal
cd ~/catkin_ws
wget http://firmware.ardupilot.org/Rover/latest/SITL_x86_64_linux_gnu/ardurover
chmod +x ./ardurover
# Copy the provided eeprom.bin in the current folder (contains the default https://github.com/ArduPilot/ardupilot/blob/master/Tools/autotest/default_params/rover.parm and the specific params from ardupilot2ros.param, note that you need to restart ardurover if you change those parameters from Mission Planner, tested with ArduRover V4.1.0-dev (833f4945)).
cp src/ardupilot2ros/other/eeprom.bin .
./ardurover --home 48.418079,-4.473487,88,0 --model rover --speedup 1

# Third terminal
cd ~/catkin_ws
mkdir -p ~/.ros/log
# Check parameters in MAVLinkDevice0.txt, warthog_cpr_sim.yaml and warthog_cpr_sim.launch if needed (all the IP addresses in these files should correspond to the computer where ardurover is running).
cp src/ardupilot2ros/other/MAVLinkDevice0.txt ~/.ros/
source devel/setup.bash
roslaunch ardupilot2ros warthog_cpr_sim.launch
# For SATURNE_simulation : roslaunch ardupilot2ros saturne.launch

Then Mission Planner can be launched and connected to 127.0.0.1:5762 (or the IP of the computer where ardurover is running). Then you can e.g. arm in Actions tab and right-click Fly to here to enter Guided mode, etc. (maybe reduce WP_SPEED and temporarily switch to Manual mode before).
