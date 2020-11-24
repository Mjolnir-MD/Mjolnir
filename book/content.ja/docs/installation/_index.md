+++
title  = "Installation"
weight = 10
+++

# Installation

Mjolnirを使用するために必要な手順を説明します。

## Prerequisities

Mjolnirを使用するために必要なものは概ね以下のようなものです。

- linux or Unix (e.g. OS X)
- C++11 compatible compiler
- git
- Make
- CMake
- Boost C++ Library

### Operating System

ファイルのパスを除けば特にlinux/unix固有の機能を使ってはいませんが、テストは
linuxとOS Xでしか行っていません。なので基本的に、Windowsでの動作は保証しません。

### C++ compiler

C++11に対応しているコンパイラが必要です。
現在は殆どの環境でデフォルトのコンパイラが対応しています。

ですが、手に入るなら、できるだけ新しいものを使用して下さい。
常にそうとは限りませんが、コンパイラが新しければより優れた最適化が行われ、
より速いプログラムになると言ってよいでしょう。

### Git

Mjolnirはtoml11という外部のライブラリを`git submodule`で管理しているため、その
ダウンロードのためにGitが必要です。

### CMake

MjolnirはCMakeを使ってビルドします。なので、CMakeが必要です。

CMakeはスクリプトの書き方の変遷が比較的速く、基本的には新しいバージョンで動く
ことを優先して開発しています。ですので、可能なら最新のバージョンを使って下さい。

### Boost

MjolnirのテストコードはBoost.Testを使用して書かれています。

Boostは、インストール済みのものが見つからなければ自動的に適切なバージョンを
ダウンロードするようにしていますが、インストール済みのものがあるならそれを使用した
方がビルドが速く済みます。

## Building

Mjolnirをコンパイルする上で必要な手順は以下の通りです。

```sh
$ git clone https://github.com/Mjolnir-MD/Mjolnir.git
$ cd Mjolnir
$ mkdir build
$ cd build
$ cmake .. # オプションは後述
$ make
$ make test
```

### Options for CMake

よく使うオプションを列挙します。

- `-DCMAKE_CXX_COMPILER=/path/to/compiler`
  - CMakeの機能です。使用するコンパイラを指定することができます。
- `-DUSE_OPENMP=(ON|OFF)`
  - デフォルトで`ON`です。
  - `ON`の場合、OpenMPが使用可能なら、OpenMPを使用したコードを含めてコンパイルします。
- `-DFIND_BOOST=ON`
  - デフォルトで`OFF`です。
  - インストール済みのBoostを探します。見つからなければ、ビルドは失敗します。
- `-DBOOST_ROOT=/path/to/boost`
  - FindBoostの機能です。Boostがインストールされているディレクトリを指定します。
- `-DSEPARATE_BUILD=(ON|OFF)`
  - デフォルトで`OFF`です。
  - `ON`の場合、分割コンパイルを行います。
  - 何度も変更してコンパイルを繰り返す開発者が時間を短縮するための機能です。
- `-DSINGLE_PRECISION_SUPPORT=(ON|OFF)`
  - デフォルトで`ON`です。
  - `OFF`の場合、`float`版のコンパイルが抑制されます。コンパイルが速くなります。
- `-DDOUBLE_PRECISION_SUPPORT=(ON|OFF)`
  - デフォルトで`ON`です。
  - `OFF`の場合、`double`版のコンパイルが抑制されます。コンパイルが速くなります。
