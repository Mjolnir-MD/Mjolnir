+++
title = "Integrators"
weight = 400
bookCollapseSection = true
+++

# Integrator

Integratorは、時間積分のアルゴリズムを意味します。

以下のアルゴリズムを選択できます。

## [BAOABLangevin]({{<relref "BAOABLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

以下の論文で提案された手法です。

- Leimkuhler B, Matthews C. Appl. Math. Res. Exp. (2013)
- Leimkuhler B, Matthews C. J. Chem. Phys. (2013)

OpenMMで推奨されているアルゴリズムです。

## [g-BAOABLangevin]({{<relref "gBAOABLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

[BAOABLangevin]({{<relref "BAOABLangevinIntegrator.md">}})と同等ですが、拘束条件を扱うことができます。

以下の論文で提案された手法です。

- Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016)

## [UnderdampedLangevin]({{<relref "UnderdampedLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

以下の論文で紹介された手法です。

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

また、以下のシミュレータで使用されている手法です。

- H. Kenzaki et al., (2011) JCTC.

## [VelocityVerlet]({{<relref "VelocityVerletIntegrator.md">}})

ニュートンの運動方程式に従い、エネルギー・体積・粒子数一定のシミュレーションを行います。
