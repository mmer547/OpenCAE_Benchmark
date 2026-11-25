# OpenFOAM v2606 / cutmesh

問題ごとに格子、初期条件、境界条件、物性、ソルバ辞書を収録しました。実行単位は `system/controlDict` のあるディレクトリです。問題名の下のディレクトリには、格子生成方式、解像度、時間刻みなどを変えた条件が入っています。96ケースと、Francis-99の形状検査入力1件を収録しています。

```bash
# OpenFOAM v2606の環境を読み込んだシェルで実行
cd sod-shock-tube/F03_100
checkMesh
if test -f system/setFieldsDict; then setFields; fi
app=$(foamDictionary system/controlDict -entry application -value)
"$app"
```

大きな `constant/polyMesh` のファイルはOpenFOAMが直接読めるgzip形式です。格子と初期条件は各ケースのメッシュに対応します。同じ問題で複数の格子を使用しているため、ケース単位で扱ってください。並列実行ではケースに含まれる `decomposeParDict` の分割数を使用し、後続ケースの初期場に使う前に `reconstructPar -latestTime` で最終時刻を復元します。

## 先行計算から初期場を作るケース

`restart.json` があるケースは先行ケースの計算終了場を使います。`source` は対象ケースから見た先行ケースの相対パスです。先行ケースを計算してから、バージョンフォルダを起点に次のように初期場を生成し、対象ケースを実行します。

```bash
python3 prepare_restart.py dfg-cylinder-re100/four_cases/DFG20dt
```

同じ格子での場の引き継ぎ、速度・圧力の倍率変更、セル中心の最近傍対応、または距離に応じた8点補間を元のケースに合わせて選びます。場の読み込みはASCII形式です。補間にはNumPy・SciPyとOpenFOAMの `postProcess` を使います。既定では先行ケースの最終時刻を読み、一部のケースでは指定時刻を使います。ダム崩壊の初期水体積分率は水柱形状から定めた入力です。

`francis-turbine-surface/` はcutmeshの形状入力検査に使ったSTLとメッシュ設定です。
