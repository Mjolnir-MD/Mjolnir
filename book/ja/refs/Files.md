# `[files]`

主にファイル入出力の設定をします。
出力フォーマットの詳細は、[Observer](Observer.md)も参照して下さい。

## 例

```toml
[files]
output.path         = "data/"
output.prefix       = "protein1"
output.format       = "xyz"
output.progress_bar = false
input.path          = "input/"
```

## 入力

`files`テーブルには、`output`と`input`の2つのテーブルがあります。

### `files.output`

`files`テーブルの中の`output`テーブルには、以下の値を設定します。

- `prefix`: 文字列型
  - ファイル名の先頭につく識別子を設定します。出力ファイルは、`{prefix}.log`のように命名されます。
- `path`: 文字列型
  - 出力ファイルのパスを設定します。そのディレクトリにファイルが出力されます。
- `format`: 文字列型
  - トラジェクトリファイルのフォーマットを指定します。
  - 粒子ごとに一つの座標情報しか格納できないフォーマットの場合、
    `{prefix}_position`と`{prefix}_velocity`の2つのファイルが出力されます。
  - `"xyz"`: [xyzフォーマット](XYZObserver.md)で出力します。
  - `"dcd"`: [dcdフォーマット](DCDObserver.md)で出力します。
  - `"trr"`: [trrフォーマット](TRRObserver.md)で出力します。
- `progress_bar`: 論理値型（デフォルトで`true`）
  - プログレスバーを表示するかどうかを選択します。コンソール出力が強制的に
    ファイルにリダイレクトされてしまう環境などで使って下さい。

### `files.input`

`files`テーブルの中の`input`テーブルには、以下の値を設定します。

- `path`: 文字列型
  - 入力ファイルを分割する際に、ルートディレクトリを設定します。
    例えば、ここに`"input"`を設定した場合、分割入力の際`./input/*.toml`が検索されます。
