+++
title  = "Glossary"
weight = 40
+++

# Glossary

Mjolnirのコードで使われる単語の用法を説明します。

## Simulator

シミュレーション全体のプロトコルです。

普通のMDシミュレーションなら指定された回数だけ時間積分を行います。

Simulated Annealingをする際は、時間ステップに応じて系の温度を変えます。

時刻に応じて力場を変更するシミュレーションでは、複数の力場を管理し、適切な
タイミングで入れ替えます。

実際の力や時間積分の計算は、Simulatorが管理する別のコンポーネントが行います。

## System

シミュレーションされる粒子の一群です。

粒子の他に、境界条件や系のパラメータ（温度など）、Topologyも管理します。

## Topology

どの粒子とどの粒子がどのような相互作用で繋がっているかを表すグラフです。

ForceFieldにおいて、結合している粒子には斥力をかけない、というような設定を
するために存在します。

## Integrator

ForceFieldに計算させた力の値を使って、時間積分を行います。

## ForceField

力場のことを指します。具体的には、PotentialとInteractionの組み合わせです。

local, global, externalの3種類があります。

### LocalForceField

あらかじめ決まった特定の粒子の間に働く力場です。

結合長ポテンシャルや結合角ポテンシャル、二面角ポテンシャルなどが代表的です。

### GlobalForceField

（電荷などのパラメータを持つ）全ての粒子の間に働く力場です。

分子間力ポテンシャル、静電ポテンシャルなどが代表的です。

### ExternalForceField

外力場との相互作用です。これのみが系全体の並進や回転を生じさせます。

空間中の特定の点へのアンカーや、系全体を覆う箱などが代表的です。

### Interaction

力場を構成する要素の一つです。ForceFieldはInteractionとPotentialから成ります。

どの粒子に力がかかるかを管理し、またPotentialを使って力の向きと大きさを計算します。

以下のような力場があるとき、

{{< katex display >}}
U(\theta) = k(\theta-\theta_0)^2,
{{< /katex >}}

実際にかかる力は

{{< katex display >}}
\nabla_i U(\theta) = \color{red}(\nabla_i \theta) \color{black} \frac{dU}{d\theta}
{{< /katex >}}

となります。この赤色の部分を実装したものがInteractionです。

### Potential

力場を構成する要素の一つです。ForceFieldはInteractionとPotentialから成ります。

ポテンシャル関数に関連するパラメータを管理し、ポテンシャル関数そのものと、その微分を計算します。

以下のような力場があるとき、

{{< katex display >}}
U(\theta) = k (\theta - \theta_0)^2,
{{< /katex >}}

実際にかかる力は

{{< katex display >}}
\nabla_i U(\theta) = (\nabla_i \theta) \color{red} \frac{dU}{d\theta}
{{< /katex >}}

となります。この赤色の部分を実装したものがPotentialです。

## Observer

エネルギーや粒子の位置、速度をファイルに出力します。

## Traits

どのような種類のシミュレーションを行うかを決めるタグ型です。

例えば、OpenMPのときにのみ必要なコードは`OpenMPSimulatorTraits`によって特殊化
されます。特殊化されなかったコードは、コンパイラが自動的に通常通りのコードを
使ってくれます。これによって実装の手間が省けています。
