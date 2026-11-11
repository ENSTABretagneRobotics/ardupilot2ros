# Copilot instructions

## Project purpose

This repository connects ROS robots to an ArduPilot Rover SITL instance.

The repository historically targets ROS 1. The current development goal is to
port it to ROS 2 Jazzy on Ubuntu 24.04 while preserving its behavior as closely
as reasonably possible.

## Repository structure

The repository contains three ROS packages:

- `ardupilot2ros`: main ArduPilot/MAVLink bridge
- `imu_conv`: IMU conversion utility
- `twist_conv`: velocity command conversion utility

Do not remove a package merely because it is inconvenient to port.

## Target environment

Target:

- Ubuntu 24.04
- ROS 2 Jazzy
- `ament_cmake`
- `colcon`

Do not attempt to retain ROS 1 Noetic/catkin compatibility unless doing so is
trivial and does not complicate the ROS 2 implementation.

## Migration principles

- Replace ROS 1 APIs with their normal ROS 2 Jazzy equivalents.
- Replace `catkin` with `ament_cmake`.
- Replace `roscpp` with `rclcpp`.
- Replace ROS 1 `tf` usage with the appropriate ROS 2 `tf2` packages where
  needed.
- Convert launch files to appropriate ROS 2 launch files.
- Convert parameter handling to ROS 2 parameters.
- Update `package.xml` files so that required packaged dependencies are
  correctly declared.
- Preserve package names, executable intent and directory organization where
  practical.
- Preserve existing functionality and algorithms.
- Do not rewrite unrelated C/C++ utility code merely for stylistic reasons.
- Preserve existing indentation and coding style wherever practical.
- Prefer small, targeted changes over broad refactoring.
- Do not silently remove functionality to make the build succeed.
- Do not disable, comment out or exclude failing targets or packages as a
  substitute for porting them.
- Do not replace working existing algorithms with substantially different
  implementations merely to simplify the migration.

## Clearpath/Warthog compatibility

Investigate the current Clearpath ROS 2 Jazzy API rather than assuming that the
old ROS 1 Warthog interfaces are unchanged.

Where the current Clearpath ROS 2 API differs from the old ROS 1 API, adapt the
launch, configuration or conversion layer where practical without unnecessarily
coupling the core ArduPilot bridge to Clearpath-specific code.

Check topic names, message types, namespaces and other relevant interfaces
against the current Clearpath documentation.

Document compatibility assumptions or integration aspects that cannot be
tested in the available environment.

Do not claim Warthog simulator compatibility for behavior that was only inferred
from documentation or compilation.

## Dependencies

You may install additional Ubuntu packages needed for investigation, building,
debugging or testing.

Prefer binary packages available through the Ubuntu 24.04 or ROS 2 Jazzy apt
repositories.

Dependencies required to build or run the ROS packages themselves must be
declared appropriately in package.xml rather than relying only on packages
manually installed in the agent environment.

Development, diagnostic and testing tools that are not dependencies of the
ROS packages do not need to be added to package.xml.

Do not modify .github/workflows/copilot-setup-steps.yml merely to make the
current session pass. If you discover an additional tool that would be useful
to preinstall for future Copilot sessions, mention it in the final report.

Use binary packages from the Ubuntu 24.04 and ROS 2 Jazzy apt repositories
wherever available.

Do not build ROS 2 or ordinary dependencies from source merely to make the
project build.

If a required dependency is unavailable as an apt package, investigate and
document the problem before introducing a source-built dependency.

Do not vendor large third-party projects merely to make the build pass unless
there is no practical packaged alternative.

If MAVLink headers or another dependency must be obtained separately, document
exactly how they are obtained and make the procedure reproducible.

A clean build should not depend on packages that happen to have been manually
installed in the development environment but are absent from the documented
installation procedure or package metadata.

## Build verification

After each logical group of changes, run a complete workspace build.

A normal verification procedure is:

    source /opt/ros/jazzy/setup.bash
    rosdep install --from-paths . --ignore-src -r -y
    colcon build --event-handlers console_direct+

Adapt the `--from-paths` path as necessary depending on the workspace layout.

If the build fails:

1. Identify the first meaningful failure rather than reacting to secondary
   errors.
2. Determine whether it is caused by source migration, ROS API migration,
   CMake, package metadata, dependency selection, generated interfaces, launch
   migration or another cause.
3. Correct the underlying problem.
4. Build again.

Continue until all repository packages build successfully or until a task
specific resource/time limit requires stopping.

Do not repeatedly run an unchanged build after a failure. A new build should
normally follow a relevant change or otherwise test a specific hypothesis.

Do not claim completion without a final successful complete-workspace build.

## Clean-build requirement

The intended final state is that a user starting with Ubuntu 24.04 and a
standard binary installation of ROS 2 Jazzy can:

1. obtain the repository;
2. install its packaged dependencies using apt/rosdep;
3. build all packages using `colcon build`.

The build must not rely on ROS 1, a previous build of this repository, or
undeclared project-specific dependencies left in the environment.

All three packages must remain part of the workspace and build successfully:

- `ardupilot2ros`
- `imu_conv`
- `twist_conv`

## Documentation

Update the repository README for ROS 2 Jazzy.

The README should contain sufficiently complete installation and build
instructions for a user starting from Ubuntu 24.04 with ROS 2 Jazzy.

It is acceptable, and useful, for the README to give explicit `apt install`
commands for the required binary dependencies.

The documented dependency installation and the dependencies declared in
`package.xml` should remain consistent.

Clearly distinguish between:

- functionality verified by compilation;
- functionality verified by execution;
- Clearpath/Warthog compatibility verified by execution;
- Clearpath/Warthog compatibility established only by inspection of interfaces
  or documentation.

Do not state that runtime or simulator functionality has been tested when only
the build has been tested.

## Scope of modifications

You may modify any tracked source code, CMake file, package metadata,
configuration file, launch file and documentation file required for the ROS 2
migration.

Changes should remain focused on the migration and compatibility objectives.
Avoid unrelated cleanup or modernization.

Do not merge the result into `main`.