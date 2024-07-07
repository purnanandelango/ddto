# Deferred Decision Trajectory Optimization

## Quadrotor motion planning examples

### `ddto-qcvx`
 - discrete-time linear system (3D double integrator)
 - convex path constraints: thrust pointing, thrust upper bound and lower bound (convexified)

### `ddto-scp`
 - continuous-time nonlinear system (3D double integrator with drag)
 - nonconvex path constraints: obstacle avoidance, speed upper bound, thrust pointing, thrust upper and lower bound