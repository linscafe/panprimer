<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

# 必須使用 hprc 第二版

**潤飾後的問題：**「我的工具必須使用 HPRC 第二版資料。」

可以；若必須以 **HPRC Data Release 2（HPRC R2）**為資料基礎，架構仍應採用「**R2 定相單倍型組裝作為 PCR 判定主體，R2 圖形泛基因組作為變異註解與座位關係層**」。R2 提供來自 232 位個體的高品質定相組裝，並規劃／提供圖形化組裝比對資源；這比 R1 更能涵蓋全球人類遺傳多樣性與複雜區域。[^1][^2]

## HPRC R2 的資料定位

HPRC R2 是預發布（pre-publication）資料釋出，包含超過 200 位個體的定序資料與高品質定相基因組；官方說明的組裝資料涵蓋共 232 位個體。組裝結合 PacBio HiFi、Oxford Nanopore 超長讀長資料及 Illumina 親子定序或 Hi-C 定相資訊，並改善組裝連續性與結構正確性。[^1]

對你的工具而言，建議將 R2 資料分成三層使用：


| 資料層 | HPRC R2 資料 | 用途 |
| :-- | :-- | :-- |
| **序列真值層** | 每位個體的 hap1 / hap2 FASTA 組裝 | 引子命中、錯配、in-silico PCR、產物序列與長度判定 |
| **座位與變異層** | 相對於 CHM13 的圖形組裝比對、VCF、assembly-to-reference alignment | 判定目標同源座位、標註 SNP、indel、SV 與引子掉擴增原因 |
| **註解與報告層** | GFF3/GTF、樣本與組裝 metadata | 基因／轉錄本座標、樣本與單倍型追蹤、統計分群與結果呈現 |

R2 的公開資源包含組裝與註解；完整的圖形化組裝比對資源則依官方資料庫持續逐步加入。官方資料頁也明確指出：R2 組裝由 Hifiasm 產生並提供對應註解。[^2]

## 參考座標選擇

HPRC R2 圖形相關 VCF 已有以 **T2T-CHM13 v2.0** 為座標基準的版本，而不是 GRCh38。[^3]

因此，建議內部使用下列原則：

- **主要座標系統：CHM13 v2.0。**
- **相容輸入：GRCh38 與 CHM13。**
- **結果輸出：同時提供 CHM13 與 GRCh38 座標，但將來源座標、lift-over 版本與不確定性明確記錄。**
- **基因與轉錄本設計：可讓使用者輸入 HGNC gene symbol、MANE transcript、CHM13 區域或 GRCh38 區域。**

這避免將 R2 的圖形變異資料強制轉回 GRCh38 時，遺失 CHM13 可表示、但 GRCh38 缺失或錯置的序列。對著絲粒附近、重複序列、大片段插入、HLA、KIR、SMN1/SMN2、CYP2D6 等區域尤其重要。

## 建議的本機資料結構

```text
hprc-r2/
├── metadata/
│   ├── samples.tsv
│   ├── assemblies.tsv
│   ├── haplotypes.tsv
│   └── data_provenance.tsv
├── assemblies/
│   ├── HG0XXXX_hap1.fa.gz
│   ├── HG0XXXX_hap2.fa.gz
│   └── ...
├── assemblies-index/
│   ├── hprc-r2-all.fa.gz
│   ├── hprc-r2-all.fa.gz.fai
│   ├── hprc-r2-all.mmi
│   └── kmer/
│       ├── k21/
│       ├── k25/
│       └── k31/
├── references/
│   ├── chm13v2.0.fa.gz
│   └── grch38.fa.gz
├── alignments/
│   ├── chm13_to_haplotype.paf.gz
│   ├── chm13_to_haplotype.chain.gz
│   └── graph/
│       ├── hprc-r2.gfa.gz
│       ├── hprc-r2.gbz
│       └── hprc-r2.dist
├── variants/
│   ├── hprc-r2.chm13.multiallelic.vcf.gz
│   └── hprc-r2.chm13.biallelic.vcf.gz
└── annotations/
    ├── chm13/
    └── assemblies/
```

這不是要求所有檔案都必須存在才可執行。第一版的必要資料是：**R2 單倍型 FASTA、組裝 metadata、CHM13 參考 FASTA，以及將 CHM13 目標區域投影至各單倍型的比對資訊**。

## R2 的索引策略

不要只建立一個將所有組裝串接的總索引。建議同時建立三種索引：


| 索引 | 建構方式 | 回答的問題 |
| :-- | :-- | :-- |
| **全域索引** | 所有 R2 haplotype contig 合併為單一 FASTA | 某個引子在整個 R2 中是否有明顯重複或非目標命中？ |
| **每單倍型索引** | 每條 hap1/hap2 組裝各自建立 k-mer 或短序列索引 | 在特定單倍型中是否多重命中或完全缺失？ |
| **座位限制索引** | 以 CHM13 目標區與兩側序列投影到每個單倍型 | 哪個命中才是預期同源座位？ |

全域索引適合快速排除重複序列；每單倍型索引則避免把「每個人都各有一份的同源座位」錯判成非唯一。座位限制索引是判斷 PCR 目標特異性所必需的層級。

## R2 版判定流程

```text
CHM13 目標區域
      │
      ├── 在 CHM13 或目標單倍型共識序列上以 Primer3 產生候選引子
      │
      ├── 對所有 HPRC R2 單倍型組裝搜尋 F / R 的精確與近似命中
      │
      ├── 使用 CHM13 ↔ haplotype 對照，標記預期同源座位
      │
      ├── 在每一個單倍型上配對 F / R 命中，建立預測 PCR 產物
      │
      ├── 使用 R2 變異圖／VCF 解釋結合位點 SNP、indel、SV
      │
      └── 依覆蓋率、唯一產物率、掉擴增與非目標擴增風險排序
```

對一組引子 $i$，建議至少計算：

$$
\text{OnTargetCoverage}_i =
\frac{\text{存在合格目標擴增子的 R2 單倍型數}}
{\text{應評估的 R2 單倍型總數}}
$$

$$
\text{UniqueProductRate}_i =
\frac{\text{僅有一個合格且目標特異產物的 R2 單倍型數}}
{\text{應評估的 R2 單倍型總數}}
$$

將每個單倍型的狀態標為 `pass`、`dropout`、`off_target`、`multi_product` 或 `uncertain`。`uncertain` 很重要：當目標區在某個組裝中斷裂、對照資訊不足，或同源投影不可信時，不應直接把它算成引子失敗。

## 引子特異性規則

以 R2 為基準時，可先採用保守的預設規則：

- 目標結合：引子全長完全匹配，或僅有遠離 3′ 端的少數錯配，且預測 Tm 合格。
- 高風險結合：引子 3′ 端最後 5 nt 出現任何錯配、引子結合區含 indel，或預測 Tm 明顯低於門檻。
- 非目標產物：正、反向引子在同一單倍型 contig 上有相向命中，且距離落在允許產物長度範圍。
- 通過條件：每個納入單倍型僅產生一個目標擴增子，且沒有可接受的非目標產物。
- 實務門檻：通用檢測可用至少 99% 的 R2 單倍型覆蓋率作為候選篩選門檻；臨床或高風險用途則應根據目標族群與實驗驗證另行訂定接受標準。

最後一點需注意：R2 中存在約兩倍於個體數的單倍型組裝，但不是所有座位都必然能在所有組裝中無歧義投影。因此，報告時應同時列出分母，例如「成功 430 / 可評估 451 個單倍型」，而非只寫成百分比。

## 圖形泛基因組的用途

若 R2 的完整圖形資料在你取得的版本中可用，建議不要讓圖直接取代 FASTA PCR 模擬；而是將它用於下列任務：

- 擷取目標基因座的**子圖**，確認引子位置是否位於泡狀變異（bubble）、插入、缺失或複雜重排。
- 找出經過某個替代分支的所有具名單倍型。
- 將掉擴增事件歸因於特定 SNP、indel 或 SV。
- 視覺化不同單倍型中產物長度差異的結構原因。
- 在高度重複或複雜基因座中，協助判斷多重命中是否來自同源等位基因、旁系同源基因或拷貝數變異。

這種分工可避免將「圖上的節點或分支數」錯當成「PCR 產物數」。PCR 是否可擴增，最終仍取決於每條實際單倍型路徑上的兩個引子結合序列、方向、距離、錯配位置與反應條件。

## R2 下載與可重現性

HPRC 指出 R2 的組裝、定序資料與註解可透過 HPRC Data Explorer、AWS S3、GitHub、AnVIL 與公共核酸資料庫取得；使用 bulk download 時，應以官方 metadata／索引檔記錄檔案位置與 accession。[^2][^1]

你的每次分析應固定紀錄：

```text
hprc_release = "R2"
download_date = "YYYY-MM-DD"
assembly_manifest_sha256 = "..."
reference_build = "CHM13v2.0"
graph_build = "..."
annotation_release = "..."
search_policy = "max_total_mm=2; max_3prime_mm=0; suffix=5"
```

由於 HPRC R2 屬於持續擴充的資料釋出，資料清單、圖形資源與註解可能更新；將 manifest 的雜湊值與工具版本納入結果，才能使每個引子設計結果可重現且可追溯。[^1][^2]
<span style="display:none">[^10][^11][^12][^13][^14][^15][^4][^5][^6][^7][^8][^9]</span>

<div align="center">⁂</div>

[^1]: https://humanpangenome.org/hprc-data-release-2/

[^2]: https://humanpangenome.org/data/

[^3]: https://zenodo.org/records/15223961

[^4]: https://hpc.nih.gov/refdb/dbview.php?id=1176

[^5]: https://www.genome.gov/sites/default/files/media/files/2023-04/hprc_publications_summary.pdf

[^6]: https://humanpangenome.org/data-use/

[^7]: https://en.wikipedia.org/wiki/Human_Pangenome_Reference

[^8]: https://projects.ensembl.org/hprc/

[^9]: https://x.com/HumanPangenome/status/1921916765577109572

[^10]: https://github.com/human-pangenomics/hpp_pangenome_resources

[^11]: https://www.genome.gov/sites/default/files/media/files/2023-04/human_pangenome_reference_consortium_executive_summary.pdf

[^12]: https://www.jinmuyu.com/download/HF02_EN.pdf

[^13]: https://github.com/Jeltje/HPRC_metadata/blob/main/hprc_metadata.tsv

[^14]: https://www.hangseng.com/cms/pws/files/pdfs/zh_hk/hab2a.pdf

[^15]: https://www.nature.com/articles/s41586-023-05896-x

