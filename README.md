# OpenCAE Benchmark

手元で検証した公開CAE問題の入力ファイルを、ソフトウェア・バージョン・問題の順に整理しています。各問題の中ではメッシュ分割、境界条件、解法などの比較条件ごとに入力を分けています。

| ソフトウェア | 主な問題 | 使用する入力 |
| --- | --- | --- |
| [OpenLB 1.9](OpenLB/1.9) | 円柱周りの流れ、自然対流、希薄気体、静止液滴、NREL停止翼など11問題 | C++、境界条件、実行スクリプト、翼のSTL・圧力採取位置 |
| [OpenFOAM v2606](OpenFOAM/v2606) | 円柱、衝撃波管、ダム崩壊、FDAノズル・ポンプ、市街地流れなど | 初期・境界条件、物性、離散化・解法辞書、cutmeshで作成した体積格子と形状 |
| [CalculiX 2.22](CalculiX/2.22) | Kirsch、NAFEMS LE1・LE10・FV52・T2～T4、Scordelis–Lo屋根 | PrePoMax 3.0.0同梱ソルバ用の `.inp` |
| [Code_Aster 16.7.9](Code_Aster/16.7.9) | 対応する構造・振動・伝熱問題と細分化条件 | Salome-Meca 2024用の `.export`、`.comm`、`.mail` |
| [Elmer 26.2](Elmer/26.2)・[音響問題](Elmer/acoustics) | 構造、振動、伝熱、音響固有モード、ダクト周波数応答 | SIF、体積・境界メッシュ、音響モデルの形状 |
| [OpenCourant 20261011](OpenCourant/20261011) | 棒の縦振動、立方体と平板の衝突、摩擦すべり、純曲げ | Starter・Engineの `.rad` と公開例からの入力作成手順 |

有限要素のケースは各条件のフォルダで実行します。CalculiXは `ccx -i 入力名`（拡張子を省略）、Code_Asterは `run_aster case.export`、Elmerは `ElmerSolver case.sif` を使います。メッシュ細分化や診断用の条件も含みます。OpenLBとOpenFOAMの準備手順は各バージョンのREADMEにまとめました。
