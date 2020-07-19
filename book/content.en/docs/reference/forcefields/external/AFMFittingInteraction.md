+++
title = "AFMFitting"
weight = 40000
+++

# AFMFittingInteraction

Flexible fitting potential to an AFM image.

It is developed in the following paper.

- T. Niina et al., JCTC (2020)

{{<katex display>}}
\begin{aligned}
U(\mathbf{r}) &= k(1 - \mathrm{c.c.}(\mathbf{r})) \\
\mathrm{c.c.} &= \frac{\sum_{p\in\mathrm{pixels}} H_p^{\mathrm{(exp)}} H_p^{\mathrm{(sim)}}(\mathbf{r})}
                     {\sqrt{\sum_{p\in\mathrm{pixels}} \left(H_p^{\mathrm{(exp)}}\right)^2}
                      \sqrt{\sum_{p\in\mathrm{pixels}} \left(H_p^{\mathrm{(sim)}}(\mathbf{r})\right)^2}} \\
H_p^{\mathrm{(sim)}}(\mathbf{r}) &= \gamma\log\left(1 + \sum_i^N \exp\left(\frac{-(x_i - x_p)^2 - (y_i - y_p)^2 }{2\sigma^2}\right)\exp\left(\frac{z_i + r_i}{\gamma}\right)\right)
\end{aligned}
{{</katex>}}

## Example

```toml
[[forcefields.external]]
interaction = "AFMFlexibleFitting"
k           = 100.0
gamma       =   1.0
pixel_x     =  10.0
pixel_y     =  10.0
length_x    =   5
length_y    =   5
sigma_x     =   2.0
sigma_y     =   2.0
z0          = 0.0
cutoff      = 5.0
margin      = 0.5
image       = [
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.5, 1.0, 0.5,
    0.0, 0.0, 1.0, 2.0, 1.0,
    0.0, 0.0, 0.5, 1.0, 0.5,
]
parameters  = [
{index = 0, radius = 1.0},
{index = 1, radius = 2.0},
{index = 4, radius = 3.0},
{index = 5, radius = 4.0},
]
```

## Input Reference

- `interaction`: String
  - Name of the interaction. Here, it is `"AFMFlexibleFitting"`.
- `k`: Floating
  - It determines the strength of the potential.
- `gamma`: Floating
  - It determines the accuracy of the softmax.
- `pixel_x, pixel_y`: Floating
  - The pixel size along each axis.
- `length_x, length_y`: Integer
  - The number of pixels along each axis.
- `sigma_x, sigma_y`: Floating
  - {{<katex>}} \sigma {{</katex>}} values in {{<katex>}} H {{</katex>}} function along each axis.
- `z0`: Floating
  - A parameter to reduce the numerical error internally. Normally, 0 is okay.
- `cutoff`: Floating
  - Cutoff length for the gaussian relative to {{<katex>}} \sigma {{</katex>}}.
- `margin`: Floating
  - Margin used in the internal neighboring list, relative to the cutoff length.
- `image`: Array of Floatings
  - The reference image. Each pixel has height in z direction.
  - The first element has (0, 0) pixel, (1, 0) pixel, ... (Lx, 0) pixel, (0, 1) pixel, ... and so on.
  - The (0, 0) pixel is the rectangular region from the origin, `(0.0, 0.0)`, to `(pixel_x, pixel_y)`.
  - The (n, m) pixel is the rectangular region from `(n*pixel_x, m*pixel_y)` to `((n+1) pixel_x, (m+1)*pixel_y)`.
- `parameters`: Array of Tables
  - `index`: Integer
    - The index of the particle.
  - `radius`: Floating
    - The radius of the particle.
