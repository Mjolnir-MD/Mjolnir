+++
title = "VelocityVerlet"
weight = 5000
+++

# VelocityVerlet

`VelocityVerlet` integrator performs constant energy simulation according to Newtonian dynamics.

## Example

```toml
[simulator]
integrator.type = "VelocityVerlet"
integrator.remove.translation = true
integrator.remove.rotation    = true
integrator.remove.rescale     = true
```

## Input reference

Some of the other parameters, such as `delta_t`, are defined in [`[simulator]`]({{<relref "/docs/reference/simulators">}}) table.

- `type`: String
  - The name of Integrator. Here, it is `"VelocityVerlet"`.
- `remove`: Table (optional)
  - `translation` and `rotation`: Boolean
    - If `true`, it removes the total translation and rotation. Otherwise, it does nothing.
  - `rescale`: Boolean
    - If `true`, it rescales all the velocities to make kinetic energy constant.
  - By default, all the fields becomes `false`.
