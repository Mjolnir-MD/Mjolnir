+++
title = "External"
weight = 3000
bookCollapseSection = true
+++

# ExternalForceField

外場との相互作用を指定します。

## [PositionRestraintInteraction]({{<relref "PositionRestraintInteraction.md">}})

粒子を空間中の点に固定します。空間中の点から一定距離の位置に固定することもできます。

## [RectangularBoxInteraction]({{<relref "rectangularbox">}})

系全体を直方体の箱の中に閉じ込めます。周期境界条件と同時には使えません。

## [ExternalDistanceInteraction]({{<relref "distance">}})

粒子と空間中の物体、壁などとの間の距離に依存したポテンシャルを指定します。

## [AFMFittingInteraction]({{<relref "AFMFittingInteraction.md">}})

AFM画像への構造のフレキシブルフィッティングを行います。

系の擬似AFM画像を生成し、参照AFM画像との相関係数に応じたポテンシャルをかけます。
