+++
title = "Integrators"
weight = 400
bookCollapseSection = true
+++

# Integrator

Integratorは、時間積分のアルゴリズムを意味します。

以下のアルゴリズムを選択できます。

## [G-JF Langevin]({{<relref "G-JFLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

以下の論文で提案された手法です。

- [Niels Grønbech-Jensen & Oded Farago, (2013) Mol.Phys. 111:8, 983-991](https://doi.org/10.1080/00268976.2012.760055)

## [BAOABLangevin]({{<relref "BAOABLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

以下の論文で提案された手法です。

- [Benedict Leimkuhler and Charles Matthews. Appl. Math. Res. Exp. (2013) 2013:1, pp. 34-56](https://doi.org/10.1093/amrx/abs010)
- [Benedict Leimkuhler and Charles Matthews. J. Chem. Phys. (2013) 138:17, 174102](https://doi.org/10.1063/1.4802990)

OpenMMで推奨されているアルゴリズムです。

## [g-BAOABLangevin]({{<relref "gBAOABLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

[BAOABLangevin]({{<relref "BAOABLangevinIntegrator.md">}})と同等ですが、拘束条件を扱うことができます。

以下の論文で提案された手法です。

- [Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016)](https://doi.org/10.1098/rspa.2016.0138)

## [GFWNPTLangevin]({{<relref "GFWNPTLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・圧力・粒子数一定のシミュレーションを行います。

以下の論文で提案された手法です。

- [Xingyu Gao, Jun Fang, and Han Wang. J. Chem. Phys. (2016) 144, 124113](https://doi.org/10.1063/1.4944909)

## [UnderdampedLangevin]({{<relref "UnderdampedLangevinIntegrator.md">}})

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

以下の論文で紹介された手法です。

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

また、以下のシミュレータで使用されている手法です。

- H. Kenzaki et al., (2011) JCTC.

## [VelocityVerlet]({{<relref "VelocityVerletIntegrator.md">}})

ニュートンの運動方程式に従い、エネルギー・体積・粒子数一定のシミュレーションを行います。
