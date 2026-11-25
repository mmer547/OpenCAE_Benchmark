# OpenLB 1.9

OpenLB 1.9r0の公式ソースを展開し、LinuxまたはWSLでC++20対応コンパイラとGNU Makeを使用します。配布アーカイブ `olb-1.9r0.tgz` のMD5は `daeb2bd2b1e1bab335aba4b73f02ce8a` です。[OpenLB 1.9.0公式配布物](https://zenodo.org/records/17899765)

```bash
export OLB_ROOT=/path/to/olb-1.9r0
cd OpenLB/1.9/dfg-cylinder-re20
bash Allrun
```

各問題の `Allrun` が必要なソースを準備し、各条件を1スレッドで順番に計算します。生成物は問題フォルダの `runs/` に入ります。

| 問題 | 条件と実行方法 |
| --- | --- |
| [定常円柱 Re = 20](dfg-cylinder-re20)・[非定常円柱 Re = 100](dfg-cylinder-re100) | 直径10・20分割で中点反射とBouzidi補間を比較。30分割はBouzidi補間 |
| [自然対流キャビティ](natural-convection-cavity) | 32・64分割、Ra = 10⁴ |
| [円管ポアズイユ流れ](pipe-poiseuille)・[同心円筒クエット流れ](annular-couette) | 中点反射とBouzidi補間の比較 |
| [周期円柱列](periodic-cylinder-array) | 固体面積率2水準、格子3水準、曲面境界2方式 |
| [移流せん断波](advected-shear-wave) | 引数 `base`、`third-order`、`product` で平衡分布を変更 |
| [希薄気体マイクロ流路](rarefied-microchannel) | 引数 `base`、`second-order`。圧力、熱クリープ、併用の3条件 |
| [円管Knudsenポンプ](knudsen-pump) | `base` は9・13分割、`pump` は17・21分割。保存補正は `conservative`（13・17）、`area`（15）、`flux`（21）。`gradient` は保存補正と解析温度勾配 |
| [静止液滴](young-laplace-droplet) | 半径15・20・25格子の圧力差と偽速度 |
| [NREL Phase VI 停止翼](nrel-phase-vi-parked) | `base` は8・6 cmの両境界条件、`fine` は5 cmのBouzidi補間、`surface` は6 cmの翼面外挿 |

停止翼のSTL、翼平面形、圧力採取位置、実測参照値は再実行用の入力です。形状を再構成するときは、報告書を `nrel.pdf`、`risoe.pdf` として `prepare.py` の隣に置き、NumPy・pdfplumberで実行します。[NREL試験報告書](https://wdh.energy.gov/docs/projects/uae6/phase_vi_test_report.pdf) ・ [Risø実験比較資料](https://backend.orbit.dtu.dk/ws/portalfiles/portal/7711119/ris_r_1465.pdf)
