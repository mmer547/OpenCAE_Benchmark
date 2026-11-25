# OpenCourant 20261011

倍精度SMP版のStarter・Engineで使用した入力です。棒の縦振動は各分割のフォルダ、他の3問題は `inputs/` にStarterとEngineのペアがあります。衝突の接触初期を短い出力間隔で追う条件は `impact/inputs-fine/` にあります。

```bash
cd axial-rod/N20
starter_linux64_gf -i rod_0000.rad -np 1 -nt 1
engine_linux64_gf -i rod_0001.rad -nt 1
```

衝突・摩擦・純曲げはAltair 2024の公開例を使用し、検証に必要な節点時刻歴、制御カードを反映しました。各問題の `input-source.json` に元の配布パスとSHA-256を記載しています。`python prepare.py friction` などで入力を再作成できます。すでに取得済みの公式ZIPを使うときは `--archive` を指定します。衝突の短時間条件は `python prepare.py impact --fine-impact` です。
