+++
title = "Files"
weight = 100
+++

# `[files]`

Here, you can control pathes and outpu file name and format.

## Example

```toml
[files]
output.path         = "data/"
output.prefix       = "protein1"
output.format       = "xyz"
output.progress_bar = false
input.path          = "input/"
```

## Input Reference

`files` table has 2 sub-tables, named `output` and `input`.

### `files.output`

- `prefix`: String
  - Filename without extension. The name output file will be `{prefix}.log` or something like that.
- `path`: String
  - The path where output files will be created. The directory must exist before running Mjolnir.
- `format`: String
  - Format of the trajectory file. One of the following.
  - `"xyz"`
  - `"dcd"`
  - `"trr"`
  - If `xyz` or `dcd` is chosen, `{prefix}_position` and `{prefix}_velocity` will be created.
- `progress_bar`: Bool (Optional. By default, `true`.)
  - If `true`, progress bar will be printed.
  - If the output is redirected to a file, Mjolnir automatically suppresses it.

### `files.input`

- `path`: String (Optional. By default, `"."`)
  - The path where Mjolnir look for files included.
