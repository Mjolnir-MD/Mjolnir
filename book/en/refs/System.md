# System

You can define positions, velocities, masses of particles, box size of a system,
and other parameters for the whole system.

It is defined as an array of tables to define several replicas of a system.
But normally, only one system is enough.

## Example

```toml
[[systems]]
attributes.temperature = 300.0
boundary_shape.lower = [  0.0,   0.0,   0.0] # lower limit of the boundary
boundary_shape.upper = [100.0, 100.0, 100.0] # upper limit of the boundary
particles = [
    {m= 1.0, pos=[ 1.000, 2.000, -1.000], vel=[ 0.100,-0.200, 0.300], name="CA", group="A"},
    # ...
]
```

## Input reference

- `boundary_shape`: Table
  - shape of the boundary. e.g. `upper` and `lower`.
  - You can omit this if the boundary type is `"Unlimited"`.
- `attributes`: Table
  - parameters of a system. e.g. reference `temperature` for a thermostat.
- `particles`: Array of Tables
  - The initial configuration of particles.

### `boundary_shape` for `"Unlimited"` boundary

You can omit the value. If it was defined, the table should be empty.

```toml
[[systems]]
boundary_shape = {} # optional.
```

### `boundary_shape` for `"Periodic"` boundary

`lower` and `upper` coordinates of the box is required.

- `lower`: Array of Floatings (length == 3)
  - The lower boundary of the box.
- `upper`: Array of Floatings (length == 3)
  - The upper boundary of the box.

```toml
boundary_shape.lower = [  0.0,   0.0,   0.0]
boundary_shape.upper = [100.0, 100.0, 100.0]
```

### `attribute`

You can define parameters of the whole system.

Some forcefields, integrators, and simulators require some specific attributes.

- `temperature`: Floating
  - The reference temperature of the system in [K].
  - NVT integrators (e.g. `BAOABLangevin`) requires this parameter.
- `ionic_strength`: Floating 
  - The ionic strength of the system in [mol/L].
  - `DebyeHuckel` potential requires this.

### `particles`

The initial configuration of the particles.

Each particle has `mass`, `position`, `velocity` (optional), `name` (optional), and `group` (optional).

```toml
particles = [
    {mass = 1.0, position = [1.0, 2.0, 3.0], velocity = [1.0, 2.0, 3.0], name = "A", group = "G"},
    # ...
]
```

Since there are many fields, you can use shortened keys.

```toml
particles = [
    {m = 1.0, pos = [1.0, 2.0, 3.0], vel = [1.0, 2.0, 3.0], name = "A", group = "G"},
    # ...
]
```

Also, you can omit velocities.
In that case, the initial velocities are generated according to Maxwell-Boltzmann distribution.

```toml
particles = [
    {m = 1.0, pos = [1.0, 2.0, 3.0], name = "A", group = "G"},
    # ...
]
```

## Defining `[[systems]]` in another file

You can split the content of a system into another file.

```toml
[files]
input.path = "./input"

[[systems]]
file_name = "system.toml"
```

```toml
# ./input/system.toml
[[systems]]
attributes.temperature = 300.0
particles = [
    {m = 1.0, pos = [1.0, 2.0, 3.0], name = "A", group = "G"},
    # ...
]
```

Of course, you can also use `include` here as another option.

```toml
[files]
input.path = "./input"

[[systems]]
attributes.temperature = 300.0
include = "initial_configuration.toml"
```

```toml
# ./input/initial_configuration.toml

# Here, we don't need [[systems]] table definition because
# the values will be expanded under `[[systems]]` table.
particles = [
    {m = 1.0, pos = [1.0, 2.0, 3.0], name = "A", group = "G"},
    # ...
]
```

## Restarting from the last snapshot of another simulation

Mjolnir saves the whole snapshot of a system in [MsgPack](https://msgpack.org/) format.

By passing a `.msg` file, you can load the whole state of the system.

```toml
[[systems]]
file_name = "restart.msg"
```
