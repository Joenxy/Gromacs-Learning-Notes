# GROMACS 全原子 MD 計算メモ

## 目的

GROMACS を用いた全原子分子動力学（MD）計算の基本的な流れを整理する。  
本メモでは、主に以下の 2 種類の体系を想定する。

1. **純水 1000 分子系**
2. **有機分子 + 水 1000 分子系**  
   例：エタノール 1 分子 + 水 1000 分子

---

## 大まかな流れ

1. MD 初期構造の作成（Packmol）
2. MD パラメーターファイルの準備
3. エネルギー最小化、平衡化計算（GROMACS）
4. 本計算（Production run, GROMACS）
5. トラジェクトリの可視化（VMD）
6. トラジェクトリの解析（GROMACS）

---

# 1. MD 初期構造の作成

## 1.1 純水 1000 分子系の場合

純水系の場合、Gaussian / RESP / GAFF2 による有機分子パラメータ化は不要である。  
必要なものは以下である。

- 水モデル構造ファイル  
  例：`SOL_GMX.pdb`
- トポロジーファイル  
  例：`topol.top`
- MD パラメーターファイル  
  例：`em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp`
- Packmol 入力ファイル  
  例：`pack.inp`

### Packmol による初期構造作成

Packmol コマンドで水 1000 分子の PDB ファイルを作成する。

```bash
packmol < input_file
```

または、入力ファイル名が `pack.inp` の場合は以下のように実行する。

```bash
packmol < pack.inp
```

実行後、水 1000 分子がランダムに配置された構造ファイルが出力される。

```text
system.pdb
```

---

## 1.2 有機分子 + 水系の場合

有機分子を含む系では、有機分子の力場パラメータが必要になる。  
既存の力場パラメータがない場合、以下のような流れで有機分子の topology を作成する。

### 有機分子パラメータ作成の流れ

1. GaussianView で初期構造を作成する
2. Gaussian により構造最適化を行う
3. 単点計算により ESP 情報を得る
4. RESP により原子電荷をフィッティングする
5. GAFF2 により結合、角度、二面角、LJ パラメータを割り当てる
6. GROMACS 用ファイルを生成する

---

## 1.3 Gaussian / GAFF2 を用いた有機分子 topology 作成

### Step 1: GaussianView による構造作成

GaussianView で分子構造を作成し、座標情報を生成する。

操作例：

```text
右クリック → Tools → Atom List
```

その後、Gaussian 入力ファイルを作成する。

```text
opt.gjf
```

---

### Step 2: Gaussian 計算の実行

以下のファイルを準備する。

```text
opt.gjf
QM.sh
sp.gjf
```

Gaussian 計算を実行する。

```bash
qsub QM.sh
```

計算後、以下のような ESP 情報を含むファイルが生成される。

```text
sp.gesp
```

---

### Step 3: GAFF2 パラメータ作成

`sp.gesp` が存在するディレクトリで以下のコマンドを実行する。

```bash
make_gaff2
```

実行後、以下の GROMACS 用ファイルが得られる。

```text
MOL_GMX.gro
MOL_GMX.top
MOL_GMX.itp
```

### 注意

`.gro` ファイルには主に以下の情報が含まれる。

- 原子座標
- 速度
- box 情報

一方、電荷、結合、角度、二面角、LJ パラメータなどの相互作用情報は主に以下のファイルに含まれる。

```text
.top
.itp
```

---

## 1.4 Packmol 用 PDB ファイルの準備

有機分子を Packmol に使うため、`MOL_GMX.gro` から `MOL_GMX.pdb` を作成する。

```bash
gmx_mpi editconf -f MOL_GMX.gro -o MOL_GMX.pdb -box 3.2 3.2 3.2
```

ここで、box サイズの単位は nm である。

---

## 1.5 Packmol による混合系の作成

Packmol 入力ファイルにより、系全体を定義する。

例えば、以下のような系を作成する。

```text
エタノール 1 分子 + 水 1000 分子
```

Packmol 入力ファイル、分子構造ファイル、水構造ファイルを同じディレクトリに置く。

```text
pack.inp
MOL_GMX.pdb
SOL_GMX.pdb
```

その後、以下を実行する。

```bash
packmol < pack.inp
```

成功すると、1001 分子が配置された初期構造ファイルが得られる。

```text
system.pdb
```

---

## 1.6 PDB ファイルを GROMACS 用 GRO ファイルに変換

Packmol で作成した `system.pdb` を GROMACS 用の `.gro` ファイルに変換する。

```bash
gmx_mpi editconf -f system.pdb -o system.gro -box X Y Z
```

ここで、`X`, `Y`, `Z` はセル長さであり、単位は nm である。

水 1000 分子程度であれば、box サイズはおよそ以下でよい。

```text
3.2 nm × 3.2 nm × 3.2 nm
```

例：

```bash
gmx_mpi editconf -f system.pdb -o system.gro -box 3.2 3.2 3.2
```

成功すると、以下のファイルが出力される。

```text
system.gro
```

---

# 2. MD パラメーターファイルの準備

GROMACS で計算を行う際には、大きく分けて以下のファイルが必要である。

## 2.1 必要なファイル

| 拡張子 | 説明 | 入力 / 出力 |
|---|---|---|
| `.pdb` | Protein Data Bank 形式の構造ファイル。MD 初期構造の作成や可視化に使用する。 | 入力 |
| `.gro` | GROMACS 固有の構造ファイル。座標、速度、box 情報を含む。 | 入力 / 出力 |
| `.mdp` | MD 計算条件を記述するファイル。温度、圧力、時間刻みなどを設定する。 | 入力 |
| `.top` | 系全体の topology ファイル。力場ファイルや分子数を指定する。 | 入力 |
| `.itp` | 分子単位または力場単位の topology 情報を含む include ファイル。 | 入力 |

---

## 2.2 `.mdp` ファイル

`.mdp` ファイルには、MD 計算条件を記述する。

例：

```text
integrator
dt
nsteps
tcoupl
pcoupl
constraints
coulombtype
rvdw
rcoulomb
```

代表的な `.mdp` ファイルは以下である。

```text
em.mdp
nvt.mdp
npt.mdp
md.mdp
```

---

## 2.3 `.top` / `.itp` ファイル

`.top` / `.itp` ファイルには、系の相互作用や分子数を定義する。

例：

```text
[ atoms ]
[ bonds ]
[ angles ]
[ dihedrals ]
[ molecules ]
```

---

# 3. エネルギー最小化、平衡化計算

Packmol などで作成した初期構造には、分子内または分子間に不自然な接触が含まれることがある。

そのまま MD 計算を始めると、力が非常に大きくなり、計算が発散する場合がある。  
そのため、通常は以下の順で計算を進める。

```text
Energy Minimization
↓
NVT Equilibration
↓
NPT Equilibration
↓
Production Run
```

---

## 3.1 エネルギー最小化

エネルギー最小化では、初期構造中の不自然な原子間接触を緩和する。

### grompp

```bash
gmx_mpi grompp -f em.mdp -c system.gro -p topol.top -o em_run.tpr
```

### mdrun

```bash
gmx_mpi mdrun -v -deffnm em_run
```

`grompp` は、以下のファイルを統合して `.tpr` ファイルを作成する。

```text
.mdp
.gro
.top
```

`mdrun` は作成された `.tpr` ファイルを読み取り、計算を実行する。

計算が完了すると、以下のファイルが得られる。

```text
em_run.gro
em_run.edr
em_run.log
em_run.trr
```

---

## 3.2 エネルギーの確認

エネルギー最小化後、以下のコマンドでエネルギーを確認する。

```bash
gmx_mpi energy -f em_run.edr -o em_potential.xvg
```

確認する代表的な項目：

```text
Potential
```

判断基準の例：

- Potential Energy が低下している
- 最大力 `Fmax` が十分小さい
- 計算が正常終了している

---

## 3.3 NVT 平衡化

NVT 平衡化では、体積を固定したまま温度を目標温度へ近づける。

例：

```text
300 K
```

### grompp

```bash
gmx_mpi grompp -f nvt.mdp -c em_run.gro -p topol.top -o nvt_run.tpr
```

### mdrun

```bash
gmx_mpi mdrun -v -deffnm nvt_run
```

温度の確認：

```bash
gmx_mpi energy -f nvt_run.edr -o nvt_temperature.xvg
```

確認する項目：

```text
Temperature
```

---

## 3.4 NPT 平衡化

NPT 平衡化では、温度と圧力を制御し、密度や box サイズを平衡化する。

例：

```text
300 K
1 bar
```

### grompp

```bash
gmx_mpi grompp -f npt.mdp -c nvt_run.gro -t nvt_run.cpt -p topol.top -o npt_run.tpr
```

### mdrun

```bash
gmx_mpi mdrun -v -deffnm npt_run
```

密度や圧力の確認：

```bash
gmx_mpi energy -f npt_run.edr -o npt_density.xvg
```

確認する代表的な項目：

```text
Density
Pressure
Box-X
Box-Y
Box-Z
Potential
```

水 1000 分子程度の小さい系では、瞬間的な圧力揺らぎは大きくなる。  
そのため、圧力は瞬間値ではなく平均値や密度の収束を確認する。

---

# 4. 本計算（Production Run）

平衡化計算が完了したら、本計算を行う。

### grompp

```bash
gmx_mpi grompp -f md.mdp -c npt_run.gro -t npt_run.cpt -p topol.top -o md_run.tpr
```

### mdrun

```bash
gmx_mpi mdrun -v -deffnm md_run
```

出力ファイルの例：

```text
md_run.xtc
md_run.trr
md_run.gro
md_run.edr
md_run.log
md_run.cpt
```

---

# 5. トラジェクトリの可視化

本計算で得られたトラジェクトリを VMD で可視化する。

```bash
vmd md_run.gro md_run.xtc
```

周期境界条件の影響で分子が切れて見える場合は、`trjconv` で補正する。

例：

```bash
gmx_mpi trjconv -s md_run.tpr -f md_run.xtc -o md_run_mol.xtc -pbc mol
```

その後、以下で可視化する。

```bash
vmd md_run.gro md_run_mol.xtc
```

---

# 6. トラジェクトリの解析

MD 計算では、各時刻における原子位置や速度の情報が得られる。  
これらを解析することで、分子集合系の構造やダイナミクスを評価できる。

ここでは、標準的な解析として以下を行う。

- 動径分布関数（RDF）
- 平均二乗変位（MSD）

---

## 6.1 RDF 解析

RDF は、ある粒子の周囲に別の粒子がどの程度存在するかを表す。

```bash
gmx_mpi rdf -f md_run.xtc -s md_run.tpr -o rdf.xvg -dt 10 -seltype mol_com << EOF
SOL
SOL
EOF
```

### コマンドの意味

```text
-f md_run.xtc
```

解析するトラジェクトリファイルを指定する。

```text
-s md_run.tpr
```

構造と topology 情報を含む `.tpr` ファイルを指定する。

```text
-o rdf.xvg
```

RDF の出力ファイル名を指定する。

```text
-dt 10
```

10 ps ごとに解析する。

```text
-seltype mol_com
```

分子の重心を用いて解析する。

---

## 6.2 Here Document について

以下の部分は Here Document と呼ばれる。

```bash
<< EOF
SOL
SOL
EOF
```

これは、通常は手入力が必要な選択項目を、シェルスクリプト内で自動入力するための記法である。

この例では、RDF の対象として水分子 `SOL` を選択している。

出力ファイル：

```text
rdf.xvg
```

---

## 6.3 MSD 解析

MSD は、分子の平均二乗変位を表す。  
拡散係数の評価に用いることができる。

```bash
gmx_mpi msd -f md_run.xtc -s md_run.tpr -o msd.xvg << EOF
SOL
EOF
```

出力ファイル：

```text
msd.xvg
```

---

## 6.4 PBC 補正後の MSD 解析

MSD 解析では、周期境界条件の影響を受ける場合がある。  
そのため、必要に応じて PBC を補正したトラジェクトリを作成する。

```bash
gmx_mpi trjconv -s md_run.tpr -f md_run.xtc -o md_run_nojump.xtc -pbc nojump
```

その後、MSD を計算する。

```bash
gmx_mpi msd -f md_run_nojump.xtc -s md_run.tpr -o msd.xvg
```

---

# 7. 拡散係数の計算

3 次元拡散では、MSD と拡散係数の関係は以下で表される。

```text
MSD(t) = 6Dt
```

したがって、MSD の傾きから拡散係数 `D` を求めることができる。

```text
D = slope / 6
```

もし MSD の単位が `nm^2`、時間の単位が `ps` の場合、

```text
1 nm^2/ps = 1.0 × 10^-6 m^2/s
```

である。

---

# 8. 注意点

## 8.1 純水系と有機分子 + 水系の違い

純水系では、Gaussian / RESP / GAFF2 によるパラメータ化は不要である。  
水モデル、例えば以下をそのまま使用できる。

```text
TIP3P
TIP4P
SPC
SPC/E
```

一方、有機分子を含む系では、その有機分子に対する topology が必要である。

既存の topology がない場合、以下のような手順で作成する。

```text
Gaussian
↓
RESP
↓
GAFF2
↓
GROMACS topology
```

---

## 8.2 `.gro` と `.top` / `.itp` の違い

`.gro` ファイルは主に構造情報を含む。

```text
座標
速度
box 情報
```

`.top` / `.itp` ファイルは主に相互作用情報を含む。

```text
原子タイプ
電荷
結合
角度
二面角
LJ パラメータ
分子数
```

---

## 8.3 GROMACS 2020 以降の cutoff scheme

古い `.mdp` ファイルでは以下のような設定が使われている場合がある。

```text
cutoff-scheme = group
```

しかし、GROMACS 2020 以降では group cutoff scheme は使用できない。  
そのため、以下の設定を使う。

```text
cutoff-scheme = Verlet
```

---

# 9. 全体フローまとめ

```text
分子構造の準備
↓
力場 / topology の準備
↓
Packmol による初期構造作成
↓
GROMACS 形式へ変換
↓
Energy Minimization
↓
NVT Equilibration
↓
NPT Equilibration
↓
Production Run
↓
Trajectory Visualization
↓
RDF / MSD Analysis
```

---

# 10. 使用ツール

| ツール | 用途 |
|---|---|
| GaussianView | 分子構造作成 |
| Gaussian | 構造最適化、ESP 計算 |
| RESP | 原子電荷フィッティング |
| GAFF2 | 有機分子力場パラメータ割り当て |
| Packmol | 初期構造作成 |
| GROMACS | MD 計算 |
| VMD | トラジェクトリ可視化 |
| Python / gnuplot | 結果の可視化 |

---

# 11. 今後追加予定

- `em.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp` の例
- TIP3P / TIP4P 水モデルの違い
- GAFF2 と OPLS-AA の違い
- RDF の見方
- MSD から拡散係数を計算する方法
- ポリマー系への応用
