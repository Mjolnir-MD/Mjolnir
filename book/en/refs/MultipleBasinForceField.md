# MultipleBasinForceField

MultipleBasinForceField connects different forcefields and makes it possible to change states between those two forcefields smoothly.

It is developed in the following paper as a method that enables conformational changes with Coarse-Grained model by employing multiple off-lattice Go potentials.

- Kei-ichi Okazaki, Nobuyasu Koga, Shoji Takada, Jose N. Onuchic, and Peter G. Wolynes PNAS (2006)

Mjolnir implements this as a generalized form, and we can use any forcefields as a "Basin" in a Multiple Basin forcefield.

In this model, the potential function $$ V_{MB} $$ is defined as a minimum eigenvalue of the matrix that appears at the left hand side of the following equation.

{% math %}
\begin{pmatrix}
V_1 + \Delta V_1 & \Delta \\
\Delta & V_1 + \Delta V_1
\end{pmatrix}
\begin{pmatrix}
c_1 \\ c_2
\end{pmatrix}
= V_{MB}
\begin{pmatrix}
c_1 \\ c_2
\end{pmatrix}
{% endmath %}

The elements of the eigenvalue, $$ c_1 $$ and $$ c_2 $$, can be interpreted as a density of the corresponding states.

As a reaction coordinate, the following value will be printed along with the energy values.

{% math %}
\chi = \log\left(\frac{c_2}{c_1}\right)
{% endmath %}

Note that this value will becomes NaN in the case of $$ 0 < \Delta $$.
To avoid this, Mjolnir makes $$\Delta$$ always negative.
$$ V_{MB} $$ depends on $$ \Delta^2 $$, so only the absolute value of $$ \Delta $$ matters.

The 3-basin case is also defined in a similar way.

{% math %}
\begin{pmatrix}
V_1 + \Delta V_1 & \Delta_{12} & \Delta_{13} \\
\Delta_{21} & V_2 + \Delta V_2 & \Delta_{23} \\
\Delta_{31} & \Delta_{32} & V_3 + \Delta V_3 \\
\end{pmatrix}
\begin{pmatrix}
c_1 \\ c_2 \\ c_3
\end{pmatrix}
= V_{MB}
\begin{pmatrix}
c_1 \\ c_2 \\ c_3
\end{pmatrix}
{% endmath %}

Here, $$\Delta_{ij} = \Delta_{ji}$$.

Currently, Mjolnir only supports 2 or 3 forcefields in MultipleBasin.

Note that all the forcefields should have the same topology.

## Example input

First, you need to define and name several forcefields to be connected.
Note that a forcefield named as `"common"` will always be applied independent from MultipleBasin forcefield.

Then, set `forcefields.type = "MultipleBasin"` in `[simulator]` table and provide parameters for MultipleBasin.

Mjolnir allows to define several MultipleBasin units.
This feature would be helpful in a case when you have several proteins that undergoes conformational changed independently from each other.

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "Unlimited"
precision     = "double"
delta_t       = 0.1
total_step    = 1000000
save_step     =   1_000
seed          = 2859805901
integrator.type = "BAOABLangevin"
integrator.gammas = [
    # ...
]
forcefields.type = "MultipleBasin"
forcefields.units = [
{basins = ["open", "close"], dVs = [0.0, -12.0], delta = 150.0},
]
# ...

[[forcefields]]
name = "open"
[[forcefields.local]]
# ...

[[forcefields]]
name = "close"
[[forcefields.local]]
# ...

[[forcefields]]
name = "common"
[[forcefields.local]]
# ...
```

## Note for AICG2+ forcefields

When AICG2+ is used with MultipleBasin, you may need to modify AICG2+ parameters slightly.
In such a case, we will have several reference structures for the same protein.
The native bond lengthes, bond angles, dihedral angles, and native contact distances could be different between reference structures and some of the native contacts could be defined only in one reference structure.

### Bond length parameter

When two native bond lengthes that correspond to the same pair of particles are too distant, the force coefficient will be moderated to decrease energy value at the transition state.

{% math %}
K'_{bond} = K_{bond} \times \mathrm{min}\left(1, \frac{E_{max}}{K_b(b_i^{(1)} - b_i^{(2)})^2}\right)
{% endmath %}

Normally, people use 100 kcal/mol as $$ E_{max} $$.

### Native contact parameter

The modification on native contacts are complicated.

There are 3 different cases.
1. The pair of particles has native contact in all the reference structures.
2. In one reference structure, the pair forms a native contact, but not in another.
3. No contact is defined in any of the reference structures.

In the case of `1.`,  the repulsive part of the native contacts will be modified.
Always the minimum native distance is used as the native distance of the repulsive part, while the attractive part is unchanged.

In the case of `2.`, the forcefield that has the native contact is unchanged.
But for the forcefield that does not have the native contact will have the repulsive part of the contact.

In the last case, `3.`, the same excluded volume interaction will be applied.
