<!-- このファイルは Markdown 記法で書かれています -->

### プログラムの概要

Lammpsのdumpファイルを解析するためのプログラム群です．
分子数などの計算条件が異なる複数のシミュレーション結果を解析する場合でも，dumpファイルが指定の書式で出力されていれば，再コンパイルは必要ありません．

### 実行手順

<Font color="aqua"> "Order" ファイルの編集 => "anal00" の実行 </font>

Orderファイル内では "\#" を行頭につけた行はコメントアウトされます．

Moleculeセクションでは分子に関する情報を行列形式で記述します．
行数（下の例では 3 ）をはじめに指定しておく必要があります．
行列は以下の手順に従って決めます．
計算系の原子を次の 2 つの条件に満たすようにグループ分けします．

1. 各グループは 1 種類の分子種しか含んでいてはいけません．
2. グループ内の原子IDは連続していなければなりません．

グループが行列の行に対応しており，1 列目がグループのタグ番号，2 列目がグループ内の分子数，3 列目がグループ内の最小の原子ID，4 列目が最大の原子IDです．
異なるグループに同一のタグ番号をつけても構いません．
また，タグ番号が昇順になるように行列をつくらなければならない，といった制約は__ありません__．
なお，タグ番号を 0 にした場合，そのグループは解析対象から除外されます．

```
$ Molecule 3 rows
0   1   1     8400
0   1   8401  16800
1   36  16801 20688
```

Atomセクションでは原子に関連した情報を行列形式で記述します．
行数はLammpsの原子タイプの数（下の例では 7 ）と同じにしてください．
2 列目以降は原子のプロパティですが，__原子の質量を必ず 2 列目に入れてください__．
それ以降はプロパティの個数を増やして自由に追加してください．

```
$ Atom 7 types with 1 property
1   12.0110
2   12.0110
3    1.0080
4   18.9980
5   15.9990
6   12.0110
7   12.0110
```

Beadセクションでは，分子内の原子とビーズの対応表ファイルの指定を行う予定です．
対応表を用いることにより，異なる分子をビーズに置換したり，ビーズへの置換方法を変えたりといったことが，プログラムを修正することなくできるようになります．

```
$ Bead
```

DumpFileセクションでは，読み込むdumpファイルの設定を行います．
読み込ませる__dumpファイルの書式__は，

<Font color="orange"> ITEM: ATOMS id type xu yu zu vx vy vz </font>

としてください．
以下，行ごとの説明です．

* Number_of_Files: dumpファイルの個数．
* Initial_Dump: 最初のdumpファイルのステップ数．
* Dump_Interval: dumpファイルのステップ間隔．
* Prefix_of_FileName: dumpファイル名のステップ数より前の部分．" " で囲むこと．
* Suffix_of_FileName: dumpファイル名のステップ数より後の部分．" " で囲むこと．
* Periodicity: 計算系の X Y Z 方向の周期境界条件の有無．

```
$ DumpFile
Number_of_Files     10
Initial_Dump        5000
Dump_Interval       5000
Prefix_of_FileName  "input_sample/atom."
Suffix_of_FileName  ".dump"
Periodicity         .true. .true. .false.
```

Analysisセクションでは，解析についての情報を記述します．
行のはじめの単語を読み込むことで，その解析が実行されます．
以降の数字は解析に必要な設定パラメータです．

```
$ Analysis
InertiaMoment 2
```

### 解析機能の追加

解析機能は追加できます．
"mod\_[解析名].f90" をつくります．
ここには，解析用の変数宣言および設定読み込み，配列割付，計算，出力サブルーチンを入れます．
"Analysis.f90" を変更します．
mod_[解析名].f90 からサブルーチンを呼び出す部分を追加します．
基本的なサブルーチンがある場合，他の解析モジュールからも使えるように "Mathematics.f90" にまとめておきます．

* InertiaMoment: 分子の慣性モーメント関連を計算し，配向秩序パラメータと回転半径を出力．

#### 構造体について

```
type cell
  double precision :: minEdge(NUM_DIMS) = 0.0d0
	double precision :: maxEdge(NUM_DIMS) = 0.0d0
	double precision :: length(NUM_DIMS) = 0.0d0
end type cell
type(cell) :: c

type molecule
  double precision :: pos(NUM_DIMS), wpos(NUM_DIMS)
	double precision :: vel(NUM_DIMS)
	double precision :: mass
  integer :: tag, minAtomID, maxAtomID
end type molecule
type(molecule), allocatable :: m(:)

type atom
  double precision :: pos(NUM_DIMS), wpos(NUM_DIMS)
	double precision :: vel(NUM_DIMS)
  integer :: typ
end type atom
type(atom), allocatable :: a(:)
```

#### 基本の機能

* Orderを読む
* DumpファイルからAtomデータを読む
*　AtomデータからMoleculeデータを計算する
*　Atomの位置データをwrapping
*　Moleculeの位置データをwrapping
* 各ステップにおけるAtomとMolecule，および計算セルのデータを構造体から参照できます．
   例） a(ID)%pos(DIM_X), m(ID)%vel(DIM_Z), c%length(DIM_Y)

#### <Font color="red">解析プログラム追加の Tutorial </font>
