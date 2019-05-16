# DCDObserver

dcd形式で位置と速度を出力します。

dcd形式は広く用いられているバイナリ形式で、多くのビューワや解析スクリプトが対応しています。

## フォーマット

dcdはバイナリフォーマットで、最初にヘッダブロックがあり、
その後各スナップショットの粒子の座標情報が続きます。

## 出力ファイル名

`[files]`の`output.prefix`で指定したファイルに、それぞれ`_position.dcd`と
`_velocity.dcd`を付け足したファイルが出力されます。

以下のような入力に対して、`data/example_position.dcd`と`data/example_velocity.dcd`
が出力されます。

```toml
[files]
output.path   = "data/"
output.prefix = "example"
output.format = "dcd"
```
