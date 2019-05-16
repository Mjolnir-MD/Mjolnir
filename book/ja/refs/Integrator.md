# Integrator

Integratorは、時間積分の手法です。

以下の手法を選択できます。

## [UnderdampedLangevin](UnderdampedLangevinIntegrator.md)

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

以下の論文で紹介されている手法です。

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

また、以下のシミュレータで推奨されている手法です。

- H. Kenzaki et al., (2011) JCTC.

## [VelocityVerlet](VelocityVerletIntegrator.md)

ニュートンの運動方程式に従い、エネルギー・体積・粒子数一定のシミュレーションを行います。
