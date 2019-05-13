## `[units]`

単位系を設定します。

`[[forcefields]]`のエネルギーのパラメータや`[[systems]]`での粒子の位置など、
ファイルから入力される数値は、この設定に従って使われます。

## 例

```toml
[units]
length = "angstrom"
energy = "kcal/mol"
```

## `units`

以下のフィールドがあります。

- `energy`: 文字列型
  - `"kcal/mol"`
  - `"kJ/mol"`
- `length`: 文字列型
  - `"nm"`
  - `"angstrom"`

質量の単位はAMU、電荷の単位は電気素量で固定されています。
現在、これらは変更することはできません。

{% hint style='working' %}
Version 1.1.0 現在、時間の単位系は指定できません。
上で指定した単位系から計算されます。
{% endhint %}

