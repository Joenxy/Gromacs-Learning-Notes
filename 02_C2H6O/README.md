Workflow：

1. MD初期構造の作成(Packmol)

2. MDパラメーターファイルの準備

3. エネルギー最小化、平衡化計算(Gromacs)

4. 本計算(Gromacs)

5.トラジェクトリの可視化(VMD)

6.トラジェクトリの解析(Gromacs)


#１．　MD初期構造の作成(Packmol)
(1-1)	Packmol用のインプットファイルを作成 
１．Gaussianview 座標情報を生成して（右クリックーTools―Atomlist）、
opt.gjfを生成する

#２．opt.gif,QM.sh,sp.gifのフィルターで以下コマンドを実行
qsub QM.sh      #sp.gespを生成する

#３．sp.gespが存在するディレクトリで以下コマンドを実行。
make_gaff2      #MOL_GMX.gro、MOL_GMX.topが得られる。groファイル、topファイルに電荷情報に座標が格納されている。


◆box作成
MOL_GMX.groファイルからbox作成に必要なMOL_GMX.pdbを作成する。
gmx_mpi editconf -f MOL_GMX.gro -o MOL_GMX.pdb -box 3.2 3.2 3.2　（㎚）

inpファイルにより系を規定して初期構造ファイル(.pdb)を作成する。下記inpファイルはエタノール1分子と水1000分子が7nm四方の立方体(box)にランダムに配置するよう指示している。
