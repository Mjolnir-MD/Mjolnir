# LocalForceField

LocalForceFieldでは、結合長や結合角、二面角など特定の粒子の間のみに働く相互作用を設定します。

## [BondLengthInteraction](BondLengthInteraction.md)

`BondLength`相互作用は、その名の通り結合長にかかる相互作用です。

## [ContactInteraction](ContactInteraction.md)

`Contact`相互作用は、2つの粒子間の距離に応じてかかる相互作用です。

なので基本的に結合長相互作用と同じ形をしていますが、結合長と違って力がかからなくなる可能性を考慮しています。
内部にポテンシャルのリストを作り、それを管理することで力のかかっていないポテンシャルをスキップすることができます。

## [BondAngleInteraction](BondAngleInteraction.md)

`BondAngle`相互作用は、3つの粒子間に形成される角度に応じてかかる相互作用です。
粒子`i`, `j`, `k`を指定すると、ベクトル`r_ji`と`r_jk`のなす角度にかかります。

## [DihedralAngleInteraction](DihedralAngleInteraction.md)

`DihedralAngle`相互作用は、4つの粒子間に定義される二面角に応じてかかる相互作用です。
粒子`i`, `j`, `k`, `l`を指定すると、`i`, `j`, `k`がなす面と`j`, `k`, `l`のなす面のなす角に対してかかります。

## [DummyInteraction](DummyInteraction.md)

`Dummy`相互作用は、何もしない相互作用です。ただし、設定したトポロジーは適用されます。

粒子`i`, `j`を指定すると、その粒子間に力をかけることなくトポロジーを設定することができます。

トポロジーの効果については、[Topology](Topology.md)を参照して下さい。
