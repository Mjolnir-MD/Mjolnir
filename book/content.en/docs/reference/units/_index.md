+++
title = "Units"
weight = 200
+++

# `[units]`

All the energy values, length, positions have the units defined here.

## Example

```toml
[units]
length = "angstrom"
energy = "kcal/mol"
```

## Input Reference

### `[units]`

- `energy`: String
  - `"kcal/mol"` or
  - `"kJ/mol"`
- `length`: String
  - `"nm"` or
  - `"angstrom"`

The unit of mass is AMU and the unit of charge is {{<katex>}} e {{</katex>}} (elementary charge).
Currently, unit of mass and charge cannot be changed.

{{<hint warning>}}
The unit of time will be calculated from the above unit system.
For example, with `"kcal/mol"` and `"angstrom"`, 1 unit-time will be approximately 49 fs.
{{</hint>}}
