<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

# 進一步說明什麼是 Exact uniqueness: k-mers

**潤飾後的問題：**「請進一步說明什麼是『以 k-mer 判定精確唯一性（exact uniqueness）』。」

**精確唯一性**是指：一段 DNA 序列在你定義的泛基因組資料庫中，是否以**完全相同的字串**只出現一次。它不允許任何鹼基錯配、插入或缺失；因此它是很快的初步篩選方法，但不能單獨證明 PCR 引子具有特異性。唯一 k-mer 通常定義為只在單一基因組位置或單一成員中出現的固定長度短序列。[^1]

## 什麼是 k-mer？

k-mer 是長度固定為 $k$ 個鹼基的連續 DNA 片段。若序列長度為 $L$，則可切出 $L-k+1$ 個重疊 k-mer。[^2]

例如，假設目標序列是：

```text
ATGCGATACG
```

若設定 $k=5$，可得到 6 個 5-mer：

```text
ATGCG
TGCGA
GCGAT
CGATA
GATAC
ATACG
```

相鄰 k-mer 只相差一個鹼基：前一個移除最左端鹼基，並在右側加入下一個鹼基。這使得整段基因組能被轉換為大量固定長度、可快速索引與計數的字串集合。[^2]

## 「精確」代表什麼？

假設你有三條單倍型組裝序列：

```text
H1: ...AACCTGATGCGATACGTTCA...
H2: ...GGGATGCGATACGAACTT...
H3: ...TTTATGCGATGCGGCCAA...
```

對 5-mer `ATGCG` 而言：


| 單倍型 | 是否有 `ATGCG` | 出現次數 |
| :-- | --: | --: |
| H1 | 是 | 1 |
| H2 | 是 | 1 |
| H3 | 是 | 1 |
| 全部合併資料庫 | 是 | 3 |

若將全部單倍型組裝直接合併為單一 FASTA，`ATGCG` 的總計數是 3，因此它**不是全域精確唯一**。但它可能是「每個單倍型中各出現一次、且都位於同源目標座位」的序列，這對設計跨單倍型 PCR 反而是好現象。

相反地，若 `GCGAT` 在 H1 的目標座位出現一次，卻也在 H1 的另一個重複區出現一次，則它在 H1 中的計數為 2。即使它在其他單倍型不存在，仍不適合作為需要高特異性的引子序列。

## Canonical k-mer

DNA 是雙股分子，因此序列與其反向互補序列在 PCR 特異性檢查中應視為同一個序列實體。

例如：

```text
序列：          ATGCG
反向互補序列：  CGCAT
```

使用 **canonical k-mer** 時，工具會在 k-mer 與其反向互補序列之間選擇固定的一種標準表示法，通常為字典排序較前者。這可以避免同一 DNA 位點因讀取方向不同而被重複計數。常見 k-mer 工具可用此作法減少冗餘。[^3]

## 在泛基因組中的層次

在泛基因組資料中，「唯一」必須明確指定比較範圍：


| 類型 | 判定方法 | 對 PCR 的意義 |
| :-- | :-- | :-- |
| **全域精確唯一** | k-mer 在所有單倍型、所有 contig 合計只出現一次 | 適合當作極高特異性的序列錨點，但通常無法做為通用人類 PCR 引子 |
| **單倍型內唯一** | k-mer 在每一條單倍型組裝中最多出現一次 | 可降低同一個人或單倍型內的非目標擴增風險 |
| **座位內唯一** | k-mer 只出現在預期基因座，可能在很多單倍型中各出現一次 | 最適合用於通用型人類 PCR 引子 |
| **群體私有唯一** | 只存在於特定單倍型、族群或等位基因，且不在其他背景出現 | 可用於基因分型、等位基因特異性 PCR 或變異檢測 |

對人類泛基因組引子設計，最實用的判定通常是：

> 候選引子或其 3′ 端 k-mer 在每個目標單倍型中，只命中一次，且命中都位於預期的同源基因座；同時，在非目標位置完全沒有命中。

這是**座位特異的精確唯一性**，而不是全域計數必須等於 1。

## 為何 k 值重要？

$k$ 越小，某個 k-mer 在基因組中重複出現的機率越高；$k$ 越大，序列辨識力通常越高，但只要出現一個 SNP 或 indel，就可能使完全匹配消失。

例如，對人類泛基因組：

- **$k=15$**：適合快速初篩，但常會命中重複序列、低複雜度區域或旁系同源基因，不適合作為唯一性證據。
- **$k=21$**：常用於定位與初步特異性篩選，在速度與辨識度間取得不錯平衡。
- **$k=25$ 至 $31$**：較適合確認長度約 20–30 nt 的 PCR 引子是否具有精確匹配的潛在非目標位點；31-mer 是常見的實務選擇之一。[^3]
- **$k$ 等於引子全長**：最直接地檢查完整引子是否精確出現；但對含常見 SNP 的結合位點過於嚴格，可能錯失可透過退化鹼基或重新設計解決的候選引子。

實務上可同時使用多個 $k$ 值，例如用 21-mer 掃描引子 3′ 端的唯一性，並用完整 20–25 nt 引子做精確命中計數。

## 與 PCR 特異性的差異

精確唯一性只回答：

> 「這個完全相同的短序列出現幾次？」

PCR 真正需要回答的是：

> 「在反應條件下，這一對引子是否可能在其他地方一起形成可擴增的產物？」

因此，以下情況即使 k-mer 顯示唯一，仍可能造成問題：

- 非目標位置與引子只有 1–2 個錯配，且錯配遠離 3′ 端。
- 正向和反向引子各自都有可結合的非目標位置，且方向相向、距離又落在允許的產物大小範圍內。
- 目標單倍型的引子結合區含 SNP 或 indel，導致某些個體無法擴增。
- 目標區位於低複雜度、重複序列、偽基因、片段重複或高度同源基因家族中。

因此，可將 k-mer 視為「快速安檢」，而將全引子比對與成對的 in-silico PCR 視為「最終核准」。

## 實作範例

若用 `meryl` 建立人類泛基因組的 25-mer 資料庫：

```bash
# 建立所有 HPRC 單倍型組裝的清單
find hprc_haplotypes -name '*.fa.gz' > assemblies.fofn

# 計算所有組裝中的 canonical 25-mer
meryl count k=25 output hprc.k25.meryl \
  $(cat assemblies.fofn)

# 取出整個資料庫中只出現一次的 25-mer
meryl equal-to 1 hprc.k25.meryl \
  output hprc.k25.global_unique.meryl
```

接著，對候選引子 `ACCTGATGCGATACGTTCAAGTGC`，可以檢查：

1. 整條 24 nt 引子在泛基因組中的精確出現次數。
2. 其 3′ 端 21-mer，例如 `TGCGATACGTTCAAGTGC`，是否在非目標位置出現。
3. 每個目標單倍型中是否皆有一次、且位於正確的同源基因座。
4. 缺少精確匹配的單倍型是否帶有引子結合位點 SNP 或 indel。

不過，上述 `meryl` 指令中的合併計數只能提供全域總數；若要判斷每個單倍型及座位的唯一性，應對每個組裝分別建置 k-mer 資料庫，或將單倍型識別資訊保留在 k-mer 索引與命中座標中。

## 對工具設計的建議

你的流程可分成三層：

1. **k-mer 過濾：** 排除低複雜度序列，以及在任一單倍型中明顯重複的候選引子或 3′ 端 k-mer。
2. **全引子比對：** 對保留的候選引子，在所有單倍型組裝中搜尋完全與近似命中，並標記錯配相對於 3′ 端的位置。
3. **引子對模擬 PCR：** 將正、反向引子的可結合位置依方向與距離配對，計算所有預測擴增子，最後依目標覆蓋率、掉擴增風險與非目標產物數量進行排序。

這種分層策略能兼顧速度與可靠性：k-mer 步驟負責快速淘汰明顯不適合的候選序列，而最終結論必須以全引子、全單倍型的 PCR 模擬結果為準。
<span style="display:none">[^10][^11][^12][^13][^14][^15][^4][^5][^6][^7][^8][^9]</span>

<div align="center">⁂</div>

[^1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC12257829/

[^2]: https://blog.csdn.net/qq_50637636/article/details/122175226

[^3]: https://blog.csdn.net/zwhszdx/article/details/148241507

[^4]: https://blog.csdn.net/ZaoJewin/article/details/140316689

[^5]: https://cloud.tencent.com/developer/article/1991230

[^6]: https://en.wikipedia.org/wiki/K-mer

[^7]: https://pmc.ncbi.nlm.nih.gov/articles/PMC7037511/

[^8]: https://escholarship.org/content/qt3p11v9rh/qt3p11v9rh.pdf?t=otc1mm

[^9]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11874746/

[^10]: https://www.airitilibrary.com/Article/Detail/U0015-0302201511580700

[^11]: https://www.bilibili.com/read/cv16360242/

[^12]: https://m.ebiotrade.com/newsf/2026-6/20260613001447302.htm

[^13]: https://cloud.tencent.com/developer/article/1613847

[^14]: https://biosciencedbc.jp/gadget/human/20160720_1_kadota.pdf

[^15]: https://pmc.ncbi.nlm.nih.gov/articles/PMC11925024/

