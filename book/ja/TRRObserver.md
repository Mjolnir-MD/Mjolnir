# TRRObserver

trr形式で位置と速度を出力します。

trr形式は広く用いられているバイナリ形式で、多くのビューワや解析スクリプトが対応しています。

## フォーマット

trrはバイナリフォーマットで、各スナップショット毎にヘッダと位置、速度、力の
情報が書かれています。

{% hint style='working' %}
基本的にtrrファイルでは長さの単位は`nm`ですが、現在のバージョンでは最初に設定した単位系になります。
{% endhint %}

## 出力ファイル名

`[files]`の`output.prefix`で指定したファイルに、`.trr`を付け足したファイルが出力されます。

以下のような入力に対して、`data/example.trr`が出力されます。

```toml
[files]
output.path   = "data/"
output.prefix = "example"
output.format = "trr"
```
