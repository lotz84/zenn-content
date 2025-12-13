---
title: "グラフからコミュニティ構造を抽出する 〜リッチフローによるグラフの時間発展〜"
emoji: "🦠"
type: "tech"
topics: ["haskell", "機械学習", "math", "graph"]
published: true
---
コミュニティ抽出とは簡単に言えばグラフにおける**ノードのクラスタリング手法**です。具体的なアルゴリズムとしてはGirvan–Newman法をはじめ様々なアルゴリズムが存在しますが、この記事では**去年（2019年）提案された新しい手法**について解説したいと思います[^1]。

[[1907.03993] Community Detection on Networks with Ricci Flow](https://arxiv.org/abs/1907.03993)

話の元になっているのはこちらの論文で、グラフを**リッチフロー**によって変形し、伸びたエッジを切断していくことでクラスタを求めるというアルゴリズムです。リッチフローという聞き慣れない言葉が出てきましたが、ちゃんと後で説明するので気にせず進めましょう。

まずは実際にグラフのクラスタリングを行う様子をアニメーションで見てみてください。

@[youtube](QlENb_XlJ_8)

アルゴリズム自体はそれほど難しくありませんが、背景を含めて理解するためには2つの理論

- **微分幾何学** （微分を使って曲がった空間を研究する分野）
- **最適輸送理論** （確率分布から確率分布への最小コストの変形（輸送）を研究する分野）

に触れておく必要があり、この記事で全てを解説するのはとても大変です。なので理論的な詳細は適宜省きながらも**アルゴリズムを理解するために必要な知識をかいつまんで解説**していきたいと思います（ちなみに線形代数の概念は説明せずに使っていくので、ベクトル空間や内積などの抽象的な定義は予め知っておくとより読みやすいかと思います）。

アルゴリズム全体がどのような考え方に基づいて作られたのかを知るためには、実は**ポアンカレ予想**について知ることが近道です。次の章では微分幾何学の概念を眺めながらポアンカレ予想がどのように解かれたのかを見ていくことにしましょう。

ポアンカレ予想
-----------
ポアンカレ予想は1904年にアンリ・ポアンカレによって提唱されてから長い間 **数学上の未解決問題**になっており、ミレニアム懸賞問題の一つとしてその解決には100万ドルの懸賞金が掛けられていました。しかし現在ではこの予想は2002年から2003年にかけてペレルマンによってarXivに投稿された3本の論文により**肯定的に解決**されています（[その後、懸賞金がどうなったかはまた別の話…](https://ja.wikipedia.org/wiki/%E3%83%9D%E3%82%A2%E3%83%B3%E3%82%AB%E3%83%AC%E4%BA%88%E6%83%B3#%E8%B3%9E%E9%87%91100%E4%B8%87%E3%83%89%E3%83%AB)）。

ポアンカレ予想の主張は以下の通りです。

> **【ポアンカレ予想】**
> コンパクトで単連結な3次元多様体は3次元球面に同相である

書かれていることを理解するだけなら実はそんなに難しくないのですが、この記事では使わない概念も多いので詳細な解説は控えておきます。ここで重要なのはこの予想が**3次元多様体**に関するものだということです。

**多様体**とは局所的にはユークリッド空間と思えるような図形（空間）のことです。例えば私達が住んでいる地球は丸いので、その表面は球面になっていますが、普段見えている景色は地面が平面であるかのように見えますよね。

![](https://storage.googleapis.com/zenn-user-upload/4gujjs1631lxvb5gmjlm9tagrgnk =500x)

このように**大域的には（全体的には）複雑な形をしているかもしれないけど、局所的に見ると（ある点の周りをすごく拡大してみると）ユークリッド空間の一部のように見える図形のことを多様体**と呼びます。局所的に見たときのユークリッド空間の次元が全て同じであることも重要で、多様体が$n$次元であるとは、局所的に見たときのユークリッド空間の次元が全て$n$だと言うことを意味します。先程の例で言えば、地球の表面は局所的にどこを見ても平面のように見えるので2次元多様体になっている、ということです（地球自体は3次元空間の中にあることに注意）。

ポアンカレ予想の証明はそれ自体を直接証明するのではなく、ポアンカレ予想を帰結として含むようなより大きな予想を解決するという形で行われました。その予想は**幾何化予想**と呼ばれています。

> **【幾何化予想】**
> コンパクトな3次元多様体は、幾何構造を持つ8つの部分多様体に分解される

幾何構造や部分多様体についても詳しくは述べませんが、大事なのはどんなコンパクトな３次元多様体でも**予め分かっている8つの基本的なパーツに分解できる**ということです[^21]。さらにそのパーツの中で"単連結"なものは3次元球面だけなので、**幾何化予想が証明されればポアンカレ予想も証明される**という関係になっています。

そしてこの幾何化予想の証明において中心的な役割を果たしたのが**リッチフロー**であり、今回のアルゴリズムの根幹となる概念なのです。

### リッチフロー
リッチフローは**曲率（空間の曲がり具合）に基づいて多様体の形を変形させていく手法**です。より正確にはリーマン多様体の計量テンソル$g_{ij}$をリッチ曲率$R_{ij}$を使って時間発展させるもので、式で書くと

$$
\frac{\partial g_{ij}}{\partial t} = -2R_{ij}
$$

のような偏微分方程式で表されます（$\frac{\partial}{\partial t}$は時間を表す変数$t$による偏微分です）。

ここでリーマン多様体・計量テンソル・リッチ曲率という３つの新しい概念が出てきたので一つずつ説明していきましょう。

まず**リーマン多様体**は微分可能多様体であって計量テンソルを備えたものです。微分可能多様体とは、多様体であって局所的に見たときのユークリッド空間上で微分ができるようなもののことを言います。微分ができると、どの方向に微分するかという情報を集めてベクトル空間を作ることができます。この空間は多様体の一点にピッタリと接した空間のように見えるので**接ベクトル空間**と呼ばれています。ここでは**微分可能多様体の各点には接ベクトル空間というベクトル空間が定義できる**ということだけ知ってもらえれば十分です。

さて微分可能多様体の各点に接ベクトル空間というベクトル空間を考えることができたので、今度はこのベクトル空間に内積の構造を入れたいと思います。微分可能多様体 $M$ の点 $p$ における接ベクトル空間 ${\rm T}_pM$ を考えて、その中の2つのベクトル ${\bm x}, {\bm y} \in {\rm T}_pM$ の内積 $\langle\cdot,\cdot\rangle$ を正定値行列 $G_p$ を用いて

$$
\langle {\bm x}, {\bm y}\rangle = {\bm x}^{\rm T}G_p{\bm y} = \sum_{i,j}g_{p,ij}x_iy_j
$$

と定義しましょう。ここで$g_{p,ij}, x_i, y_j$はそれぞれ$G_p, {\bm x}, {\bm y}$の成分を表しています。接ベクトル空間は多様体の各点に存在するので、$g_{p,ij}$も各点上に定義されていますが、$g_{p,ij}$が$p$を変数とする関数として滑らかに変化するとき、それ全体を$g_{ij}$で表し、**計量テンソル**と呼びます。内積の構造は接ベクトル空間におけるベクトルの長さや向きといった量を決定づけるので、計量テンソルは（微分可能多様体全域で整合的に定義された）ベクトルの長さや向きに関する構造だと言えます。本当にざっくりした表現ですが**計量テンソルは多様体の形に関する量**だと考えられるでしょう。

いよいよ**リッチ曲率**について見ていきましょう。

曲率は空間の曲がり具合を測る量です。空間が曲がっているとはどういう状態でしょうか？例えば平面上にあるベクトルを平行移動することを考えましょう。この時、当たり前ですがどのように平行移動したとしてもベクトルの向きはずっと一定です。しかし今度は球面上でベクトルを平行移動することを考えましょう。下図のように球面上にA,B,Cの3点を取り、接平面上のベクトルをAからCに平行移動するパターンと、AからBを経由してCに平行移動するパターンを考えると両者の結果は一致しません。このように**曲がった空間の上では平行移動した結果が経路に依存してしまう**のです。

![](https://storage.googleapis.com/zenn-user-upload/1b9l08cyeg1486rgoc5qpa9wfger =400x)

この現象を逆手に取って、どの方向に進めばどれだけベクトルがずれるのかを測定することで空間の曲がり具合を調べることができます。この時に現れるのが**断面曲率**と呼ばれる曲率です。断面曲率の正確な定義を述べることは避けますが、断面曲率は以下のような場面で現れます。

ある接ベクトル空間 ${\rm T}_xM$ の二つの単位ベクトル ${\bm v}, {\bm w}$ と十分小さい正の実数 $\delta, \varepsilon \in \mathbb{R}$ を使って**辺の長さが$\delta, \varepsilon$の平行四辺形を作る**ことを考えましょう。まず$x$から$\varepsilon{\bm w}$だけ移動した点を$x'$とします。次に$x$から$\delta {\bm v}$移動した先にある点を$y$とおいて、さらに$x$から$\delta{\bm v}$に沿って${\bm w}$を平行移動したベクトル${\bm w'}$を使って$y$から$\varepsilon{\bm w'}$移動した点$y'$を考えます[^2]。

![](https://storage.googleapis.com/zenn-user-upload/ep85q7su60plrus1xw2vk1ek63as =600x)

（左: 平坦な空間での平行四辺形、右: 曲がった空間で同様の手順を行ったもの）

この手順により、曲がっていない平坦な空間では綺麗な平行四辺形が作れます。しかし曲がった空間ではそうなるとは限らず、実際$xy$の距離$d(x, y)$と$x'y'$の距離$d(x',y')$が一般的には一致しません[^3]。$x'y'$の距離は

$$
d(x', y') = d(x, y)\left(1-\frac{\varepsilon^2}{2}K({\bm v}, {\bm w})+O(\varepsilon^3+\varepsilon^2\delta)\right)
$$

と書けることが知られており[^4]、この式に現れる$K$が断面曲率です。ちなみに$O$はランダウの記号であり、高次の項は$\delta, \varepsilon$が十分小さい時は無視できることが分かります。

断面曲率はどの方向${\bm v}$にどんなベクトル$w$を平行移動するのかという２つのベクトルを指定する必要がありましたが、1つのベクトル${\bm v}$だけを使ってその方向での**平均的なベクトルの変化率によって空間の曲がり方を調べる**こともできます。

ある接ベクトル空間 ${\rm T}_xM$ の単位ベクトル ${\bm v}$ と十分小さい正の実数 $\delta, \varepsilon \in \mathbb{R}$ を考え、点$x$を中心とする半径$\varepsilon$の球体を$\delta{\bm v}$に沿って平行移動させることを考えます。

![](https://storage.googleapis.com/zenn-user-upload/gixgre1uxiz65logcxr3lqgtyquh =300x)

この時、球体の点が移動する距離の平均は

$$
d(x, y)\left(1-\frac{\varepsilon^2}{2N}R({\bm v}, {\bm v})+O(\varepsilon^3+\varepsilon^2\delta)\right)
$$

と書けることが知られており[^4]、ここで現れる$R$がリッチ曲率と呼ばれるものです。つまり**リッチ曲率はある方向${\bm v}$に進んだ時の空間の大きさの変化に関する量**だと言えるでしょう。リッチ曲率は二階テンソルになっておりその成分表示を$R_{ij}$と書くことにします。

これでリッチフローの方程式を構成する概念が出そろいました。再掲すると

$$
\frac{\partial g_{ij}}{\partial t} = -2R_{ij}
$$

という式でした。実はリッチ曲率は適切な座標系を考えると

$$
R_{ij} = -\frac{1}{2}\Delta g_{ij} + {\rm lower\ order\ term}
$$

と書けることが知られており[^6]、組み合わせると

$$
\frac{\partial g_{ij}}{\partial t} = \Delta g_{ij}
$$

だいたいこのような式になります（lower order termは無視しました）。これは熱伝導方程式に似た形になっていますね。熱伝導方程式は熱が伝わっていく現象を記述する方程式で、例えば長い金属の棒の中心を熱した後に放置することを考えると、熱はどんどん周りに伝わっていってなだらかな分布になっていくと思います。

リッチフローは熱の代わりに計量テンソルに関する熱伝導方程式のように思えるので、**リッチフローは複雑な計量テンソル、つまり多様体の形をより滑らかなものに変形していく過程**だということが分かります。

ここからはリッチフローの具体例を見ていきましょう。

まずは**2次元多様体のリッチフロー**の例を見てみましょう。以下の図は上から下に時間が流れており、円柱によって連結された３つの球面がどんどん単純になり単純な球面に収束していく様子が分かります。

![](https://storage.googleapis.com/zenn-user-upload/wyvjdd22pel5gfmsto5tnhfrrs69 =200x)

（画像出典: [Ricci flow - Wikipedia](https://en.wikipedia.org/wiki/Ricci_flow)）

実はコンパクトな2次元多様体はリッチフローによる変形で曲率が一定に収束することが知られています。曲率が正の時は楕円幾何、曲率が0の時はユークリッド幾何、曲率が負の時は双極幾何と呼ばれていて、2次元多様体の幾何構造はこの3種類に分類されるのです[^7]。

次にコンパクトな**3次元多様体のリッチフロー**の例を見てみましょう。2次元多様体の場合は3種類のどれかに収束しましたが、3次元多様体の場合はそれほど単純ではありません。例えば以下のネック・ピンチと呼ばれる構造を考えます。

![](https://storage.googleapis.com/zenn-user-upload/haa19b0esmtjbzbqmyg82ln3qmax =400x)

これは一見すると先程の2次元多様体の時と同じ様な図ですが、今度は3次元多様体です。これは3次元球面と3次元球面を（円柱ではなく）2次元球面柱（2次元球面と単位区間の直積）で滑らかに繋いだような多様体をリッチフローによって時間発展させたもの（を平面に投影した絵）になっています。実は2次元多様体の場合、円柱の曲率は0なので球面同士をつなげている間の構造自体に目立った変化はなかったのですが、3次元多様体の場合2次元球面の曲率が正なので間の構造は**どんどん細くなっていき1点に潰れてしまう現象**が起こるのです。

幾何化予想はコンパクトな3次元多様体をリッチフローによって変形することで8つの幾何構造に分解できることを示すというアプローチで解決されました。リッチフローによる変形の際に現れる上述のような特異点をペレルマンは手術と呼ばれる方法で取り除くことを行いました。手術は細く長くなった部分が一点に潰れる前に切断し、分かれた二つの断面をそれぞれ滑らかに塞ぐという手順で行われます。この操作によって**多様体がどんどんより単純な部品に分解されていく**のです。

![](https://storage.googleapis.com/zenn-user-upload/dppg75des96u0m03thane2q9cm7u =400x)

さてポアンカレ予想の話はこれぐらいにして、そろそろコミュニティ抽出アルゴリズムの話に戻りましょう。

:::message
微分幾何学に関連する概念の正確な定義やリッチフローの話を更に知りたい人は以下の文献が読みやすかったので参考までに載せておきます。

* [Nick Sheridan - Hamilton’s Ricci Flow](https://homepages.warwick.ac.uk/~maseq/topping_RF_mar06.pdf)
* [Peter Topping - Lectures On The Ricci Flow](https://homepages.warwick.ac.uk/~maseq/topping_RF_mar06.pdf)
:::

Ollivier-Ricci曲率
------------------
幾何化予想解決のアイデアはリッチフローによって3次元多様体を変形し基本的な部品に分解するというものでした。この考え方をグラフにそのまま適応してみましょう。つまり**グラフをリッチフローによって変形することでコミュニティという基本的な部品に分解していく**というアイデアです。そのためには微分可能多様体の上に定義されていたリッチ曲率という概念をグラフに対しても定義する必要が出てきます。しかしいきなりグラフ上でリッチ曲率を考えるのは難しいので**まずは距離空間を経由する**ことにします[^13]。

リッチ曲率はベクトル場や共変微分など微分幾何学の言葉で定義された概念です。しかし上述のリッチ曲率が出てきた例においては2点間の距離とそれぞれの点を中心とする球体の平均移動距離によって特徴づけられていました。Ollivierは**この式を逆に解くことで確率分布を備えた距離空間におけるリッチ曲率を定義する**ことを考えたのです[^8][^9]。

$(X, d)$ を距離空間とし $X$ の各点 $x$ に添字付けられた確率分布の族 $\left(m_x\right)_{x\in X}$ が与えられているとします（ $m_x$ は $x$ を中心とするガウス分布のようなものを考えておくと良いでしょう）。この時 $X$ 上の2点 $x, y$ 間の **Ollivier-Ricci曲率**$\kappa$を

$$
\kappa(x, y) = 1 - \frac{W_1(m_x, m_y)}{d(x, y)}
$$

と定義します。ここで$W_1$はWasserstein距離と呼ばれるもので、2つの確率分布の間の距離を測るものです。この式の意味はまた戻ってきて確認するとして、まずはWasserstein距離について知るため手短に最適輸送理論の話をしていきたいと思います。

### 最適輸送理論
最適輸送理論の発端となる問題は18世紀後半にフランスの数学者Mongeによって提起されました。それはある砂の山を別の地点まで運んで再び砂の山を作る時、どのように運べば輸送コストを最も小さくできるのかという問題でした。これは砂の山を確率分布だと思えば、**ある確率分布から別の確率分布への変形を考えた時にどの様な変形を考えればコストを最小に押さえられるかという問題**と解釈することができます。そして実はこの時の最小コストが2つの確率分布の間の距離になるのです。

![](https://storage.googleapis.com/zenn-user-upload/r8hdfh4350jd5yf2voi1lvyfajf8 =400x)

最適輸送理論は一般的には測度の言葉で記述されるものですが、この記事では有限集合上の確率分布しか扱わないので確率ベクトル（値が全て正の実数であり合計が1になるようなベクトル）を使って説明していくことにします。

Mongeは前述の問題を数学的に定式化しましたが、Mongeの方法では解が存在しない可能性があったり、有限集合では解の探索空間が組み合わせ的に大きくなってしまう等とあまり性質がよくありません。これに対しKantorovichはMongeの問題の制約を幾分か緩めた以下のような定式化を考えました。

二つの有限な点集合$X=\{x_1, \dots, x_m\}, Y = \{y_1, \dots , y_n\}$とそれぞれの集合上の確率分布に対応する確率ベクトル ${\bm a} = (a_1, \dots, a_m), {\bm b} = (b_1, \dots, b_n)$、そして点$x_i$から点$y_j$への輸送コストを表す関数$c(x_i, y_j)$ を考えます。この時、以下のような0以上の実数を値に持つ$n\times m$行列の集合

$$
U({\bm a}, {\bm b}) = \{P\in\mathbb{R}^{n\times m}_{\geq 0} | P\mathbb{1}_m = {\bm a}, P^{\rm T}\mathbb{1}_n = {\bm b}\}
$$

の中から、

$$
\min_{P \in U({\bm a}, {\bm b})} \sum_{i,j}c(x_i, y_j)p_{ij}
$$

を満たす$P$を見つけるという問題です。ここで$1_n$は全成分が$1$の$n$次元ベクトルを表し、$p_{ij}$は行列$P$の$i,j$成分を表しています。

$P$は行方向に足し合わせると${\bm a}$が、列方向に足し合わせると${\bm b}$が出てくる行列になっていて、${\bm a}$と${\bm b}$のカップリング行列と呼ばれます。$P$の全要素を足し合わせると1になることを考えると、$p_{ij}$の値は地点$x_i$から$y_j$にどれぐらいの量の"砂を運ぶか"と解釈することができ、Kantorovichの問題は **${\bm a}$から${\bm b}$への総輸送コストを最小にする$P$を見つける問題** だと考えることができます。この問題は線形の等式・不等式制約の下、線形な式の最小化を考えているため**線形計画法によって解くことが可能**です。

![](https://storage.googleapis.com/zenn-user-upload/e8892o9rqw3nc0ycdmoq26d28zcz =400x)

Kantorovichの問題において$n=m$すなわち${\bm a}$と${\bm b}$が同じ次元であり、$c(x_i, y_j)$がある距離関数$d$の$q$乗で表される時（ただし$q\geq 1$）、

$$
W_q({\bm a}, {\bm b}) = \left(\min_{P \in U({\bm a}, {\bm b})} \sum_{i,j}d(x_i, y_j)^qp_{ij}\right)^{1/q}
$$

と置くと$W_q({\bm a}, {\bm b})$は$n$次元の確率ベクトル空間の距離になることが知られています。これは元の集合$X$や$Y$に属する点の間の距離ではなく、その上に定義される**確率分布の間の距離**です。この距離は**Wasserstein距離**と呼ばれています。

Kantorovichの問題は線形計画法で解けるためかなり効率的に解ける問題ではありますが、実はこの後グラフのエッジ数に比例する大量のWasserstein距離を計算する必要があることを考えると、さらなる計算の高速化を準備しておきたいところです[^10]。そのためKantorovichの問題にさらにエントロピーの制約による正則化を加えることで、**Shinkhornアルゴリズムとよばれる手法によってより高速に問題が解けるようになる**ことを使いたいと思います。ただし高速に解けると言っても元の問題に正則化を加えたものなので答えが正確に一致するとは限らない近似的な解法です。

カップリング行列$P$のエントロピー$H(P)$を

$$
H(P) = -\sum_{i,j}p_{ij}(\log p_{ij}-1)
$$

と定義します。考えるのはこのエントロピーによってKantorovichの問題を正則化した

$$
\min_{P \in U({\bm a}, {\bm b})} \left(\sum_{i,j}c(x_i, y_j)p_{ij} - \varepsilon H(P)\right)
$$

という問題です。エントロピー$H(P)$がカップリング行列$P$がどれだけ乱雑か/分散してるかを表す指標だと考えると、**正則化付きのKantorovichの問題は$\varepsilon$程度カップリング行列がぼやけることを許容した問題**だと考えることができるでしょう。実際以下の図は$\varepsilon$を変化させた時の解となるカップリング行列の変化を表していますが、**$\varepsilon$が小さくなればなるほどよりシャープな解に近づいている**ことが分かります。同時に$\varepsilon$が小さくなるほどKantorovichの問題の解に近づきます。

![](https://storage.googleapis.com/zenn-user-upload/9s74zgofrju0lbcq7hpjqd5xhy3l =500x)
（画像出典: [[1803.00567] Computational Optimal Transport](https://arxiv.org/abs/1803.00567)）

正則化付きのKantorovichの問題を考える大きな理由は、**この問題を解く高速なアルゴリズムが知られている**ことです。それは**Shinkhornアルゴリズム**と呼ばれています。

まず正則化付きのKantorovichの問題の解となるカップリング行列$P$は

$$
P = {\rm diag}({\bm u})K{\rm diag}({\bm v})
$$

と表せることが知られています。ここで$u, v$は非負の実数を要素に持つベクトルで${\rm diag}$は与えられたベクトルを対角成分に持ち残りの成分が0であるような行列を返す関数です。また行列$K$は$ij$成分が

$$
K_{ij} = \exp\left(-\frac{c(x_i, y_j)}{\varepsilon}\right)
$$

であるようなものとします。

つまり解となるカップリング行列$P$を見つけるためには行列$K$をスケーリングする二つのベクトル$u, v$を見つければ良いわけです。そしてこれは以下の計算を反復的に行うことで求まります。

$$
{\bm u}^{n+1} = \frac{\bm a}{K{\bm v}^n},\ \ {\bm v}^{n+1} = \frac{\bm b}{K^{\rm T}{\bm u}^{n+1}}
$$

ただしベクトルの割り算は成分毎の割り算を行うことを表しています。あとはこの計算を適当な初期値から初めて収束するまで繰り返せば、目当ての${\bm u}, {\bm v}$を計算することができます。

ここまで最適輸送理論について一つの証明もなく必要な事実を淡々と述べてきましたが、ここでの話は [Computational Optimal Transport](https://arxiv.org/abs/1803.00567) を元にしているので詳細など詳しく知りたい方はぜひこちらの文献を参照してください。

Wasserstein距離の定義とその計算方法が分かったので、そろそろリッチ曲率の話に戻りましょう。

### グラフ上のリッチ曲率
もう一度Ollivier-Ricci曲率の定義を思い出しましょう。

$$
\kappa(x, y) = 1 - \frac{W_1(m_x, m_y)}{d(x, y)}
$$

$W_1$が$p=1$のWasserstein距離であることはもう分かりますね。$m_x$として点$x$を中心とするガウス分布のような確率分布を考えるとWasserstein距離$W_1(m_x, m_y)$は点$x$周辺の空間を点$y$周辺の空間に移動する時のコストを表していること解釈することができます。Ollivier-Ricci曲率はこの値と$x,y$の距離$d(x,y)$を比較することで**移動距離に対して空間が大きくなっているか小さくなっているかを測ろうとしている**のです。実際、距離空間としてリーマン多様体を持ってきて$m_x$を点$x$を中心とする$\varepsilon$球の部分にだけ一様な値を取る確率分布を考えると、Ollivier-Ricci曲率は$\varepsilon, d(x, y)$を0に近づけることでRicci曲率に収束することが分かります。

Ollivier-Ricci曲率はその定義が明示的にも暗黙的にも何らかの極限に依存していないので、離散的な距離空間であってもOllivier-Ricci曲率を考えることができます。つまり**グラフに対してもOllivier-Ricci曲率を定義することができる**のです（ここからはグラフを距離空間とみなせるように正の重みがついたループも多重辺もない無向グラフを考えます）。

グラフのノード$x$に対するグラフ上の確率分布$m^\alpha_x$を以下のように定義しましょう。

$$
m^\alpha_y =
\begin{cases}
\alpha & y = x \\
\frac{1-\alpha}{C}\exp\left(-d(x,y)\right) & y \in \pi(x) \\
0 & {\rm otherwise}
\end{cases}
$$

ここで$\alpha \in [0, 1]$は確率分布のパラメータ、$\pi(x)$はノード$x$とエッジを共有しているノードの集合、$C$は合計を1にするための規格化定数で$C = \sum_{x\in\pi(x)}\exp\left(-d(x,y)\right)$です。定義からこれは**ノード$x$の上に確率$\alpha$が乗っていて、残りの確率は距離に応じて減衰するように隣接ノードに配られているような確率分布になっている**ことが分かるでしょう。

さて、グラフ上の確率分布が定義できたのでこの確率分布を使って二つのノード$x,y$間におけるOllivier-Ricci曲率を定義しましょう。

$$
\kappa^\alpha(x, y) = 1 - \frac{W_1(m^\alpha_x, m^\alpha_y)}{d(x, y)}
$$

式自体は変わりませんが、リッチ曲率が確率分布のパラメータに依存していることを明示的に表すようにしました。$d(x, y)$はノード$x, y$間の最短距離であり、Wasserstein距離を計算する際の輸送コストもグラフ上の最短距離を用いることとします。

実は$\alpha$に対して以下のような極限を取ると、グラフの直積における曲率が元のグラフの曲率から計算できるなど良い性質を持つことが知られています[^14]。

$$
\kappa(x, y) = \lim_{\alpha \rightarrow 1}\frac{\kappa^\alpha(x, y)}{1-\alpha}
$$

実装する際に極限を計算するのは難しいですが、$\kappa^\alpha(x, y)$は$\alpha$について高々3つの区間からなる区分線形であることが知られている[^15]ので十分大きな$\alpha$を使えば問題ありません。グラフ上のリッチ曲率としてはこの$\kappa(x, y)$を採用することにします。

ここからは実際の例でリッチ曲率がどうなるか見てみましょう。

まずは以下のような立方体状のグラフを考えましょう。エッジの長さは全て1とします。

![](https://storage.googleapis.com/zenn-user-upload/xzuvw6zqsrxoojreejofg62tm8jj =200x)

このグラフ上の確率分布$m_A, m_B$は以下のように考えることができます。

![](https://storage.googleapis.com/zenn-user-upload/qllvbbgh0knhnrt9ygglk0czhbfi =500x)

 $\alpha$と$\frac{1-\alpha}{3}$の大小によって結果は変わってきますが、$\alpha > \frac{1-\alpha}{3}$と仮定すると、Wasserstein距離$W_1(m_A, m_B)$は

$$
\begin{matrix}
W_1(m_A, m_B) &=& \frac{1-\alpha}{3} + \frac{1-\alpha}{3} + \left(\alpha - \frac{1-\alpha}{3}\right) \\
&=& \frac{1-\alpha}{3} + \alpha \\
&=& \frac{1+2\alpha}{3} \\
\end{matrix}
$$

となります。ここで輸送コストはグラフ上のノード間の最短距離です。Ollivier-Ricci曲率は

$$
\begin{matrix}
\kappa^\alpha(A, B) &=& 1 - \frac{W_1(m_A, m_B)}{d(A, B)} \\
&=& 1 - \frac{1+2\alpha}{3} \\
&=& \frac{2}{3}(1-\alpha) \\
\end{matrix}
$$

と計算できます。更に$1-\alpha$で割れば（極限を取る必要もなく）、

$$
\kappa(A, B) = \frac{2}{3}
$$

となることが分かります（ちなみに同様の議論で$\alpha \leq \frac{1-\alpha}{3}$の場合は、$\kappa^\alpha = 2\alpha$となります）。

さらに立方体の対称性より全てのエッジにおけるリッチ曲率が同じ値になります。**立方体が３次元球を近似していると考えると曲率が正になるのは直感に合いますね**。

次に以下のようなグラフを考えてみましょう。再びエッジの長さは全て1とします。

![](https://storage.googleapis.com/zenn-user-upload/7m8pgzgv9fje0254vgpnlcwnkjk5 =400x)

A, Bは同じコミュニティに属するようなノードで、C, Dはコミュニティ間を接続するようなノードになっていますね。

前の例と同様の計算でA, B間のリッチ曲率は

$$
\kappa(A, B) = \frac{3}{2}
$$

となり正の値になることが分かります。そしてC, D間のリッチ曲率は

$$
\kappa(C, D) = -2
$$

となり負の値になることが分かります。

一般にコミュニティを形成しているノード間のリッチ曲率はノード間の繋がりが密なので周辺にショートカットに使える道が多くWasserstein距離の方がノード間の距離よりも小さくなる傾向があります。つまり**コミュニティを形成しているノード間のリッチ曲率は正になる傾向があります**。

反対にコミュニティ間を繋いでいるブリッジのようなノードのリッチ曲率は、周辺ノードの距離はそのブリッジを経由することが多くなりWasserstein距離の方がノード間の距離よりも大きくなる傾向があります。つまり**コミュニティ間のブリッジのノードのリッチ曲率は負になる傾向があります**。

リッチフローはリッチ曲率にマイナス符号をつけた向きに空間を変形させていくのでした。つまりグラフ上のリッチフローを考えると**同じコミュニティに属するノード間の距離はどんどん小さくなり、コミュニティの間を繋ぐようなエッジはどんどん長くなってコミュニティを引き離していく**ような挙動が期待できそうです。

それでは、いよいよグラフ上のリッチフローを見ていきましょう。

グラフ上のリッチフロー
-----------------
論文[^16]によればグラフのリッチフローは

$$
\frac{dw_{xy}}{dt} = -\kappa_{xy}w_{xy}
$$

という微分方程式を離散化することで得られます。

$$
w^{n+1}_{xy} = w^n_{xy} - \varepsilon\kappa^n_{xy}w^n_{xy}
$$

ここで$w^n_{xy}$は$n$ステップ目のノード$xy$間のエッジの重みであり、$\kappa^n_{xy}$は$n$ステップ目のノード$xy$間のリッチ曲率です。もし$\varepsilon$が1であれば右辺は$m_x, m_y$間のWasserstein距離に等しくなります。

さっそくリッチフローを使ったグラフの時間発展を考えても良いのですが、このままだとコミュニティ構造がどんどん小さくなってしまうので、正規化することを考えましょう。離散化したリッチフローの各ステップでの計算の直前にノード$xy$間のエッジの重み$w_{xy}$を以下のように更新します。

$$
w^n_{xy} \leftarrow d^n(x,y)\frac{|E|}{\sum_{xy}d^n(x,y)}
$$

ここで$d^n(x,y)$は$n$ステップ目のノード$xy$の最短距離、$|E|$はグラフの全エッジの個数、分母の和は全てのエッジについて取ることとします。つまり全エッジの重みの和が一定になるように維持するのです[^17]。

最後にコミュニティを分割するためのエッジの切断（手術）を考えましょう。といってもやることは単純で、元の論文ではリッチフローのステップを5回計算する毎に上位5%の重みのエッジを切断しているようです。しかし同じ論文に通常は10~15回に1回行うという記述もあり、切断の条件は実験しながらチューニングしていくのが良さそうです。GitHubにある筆者らの実装ではModularityというグラフの指標に基づいて良い感じにエッジの切断を行うような実装がされてたりもします[^19]。

ここまでの手順をまとめると以下のようなアルゴリズムになります。

**【離散化したリッチフローのアルゴリズム】**
重み付きのグラフ$G$と収束判定のしきい値$\delta$が与えられる

1. グラフの重みを正規化する
2. グラフのリッチフローに基づいてグラフの重みを更新する
3. 条件を満たしているエッジを切断する
4. 全てのリッチ曲率の変化率が$|\kappa^{n+1}_{xy} - \kappa^n_{xy}| < \delta$を満たすまで1-3を繰り返す

最終的にエッジの重みがリッチフローによって更新されたグラフが出力されます。

Haskellによる実装
---------------
ここまで散々数学の話をしてきましたがZennはエンジニアのための情報共有コミュニティなので、グラフ上のリッチフローを実際に実装しないことには終われません。ここからはHaskellを使ってグラフ上のリッチフローを実装していくことにしましょう[^19]。

グラフ上のリッチフローを実装するために必要な計算は

- ノード間の最短距離の計算
- グラフ上のWasserstein距離の計算

の二つです。実は2日前に公開した記事でノード間の最短距離を計算する方法を紹介しているので、まずはこの記事と同じように`Semiring`と`Vector n a`で線形代数を実装していきましょう。

https://zenn.dev/lotz/articles/9c9ca0708b035b

### 線形代数

まずはベクトルに関する計算を行う関数の実装です。

```hs
import Data.Vector.Sized (Vector)
import qualified Data.Vector.Sized as V

-- | ベクトルのスカラー倍
(*^) :: Num a => a -> Vector n a -> Vector n a
a *^ v = fmap (a*) v

-- | ベクトル間の足し算
(^+^) :: Num a => Vector n a -> Vector n a -> Vector n a
(^+^) = V.zipWith (+)

-- | ベクトル間の引き算a
(^-^) :: Num a => Vector n a -> Vector n a -> Vector n a
(^-^) = V.zipWith (-)

-- | ベクトルのL2ノルム
norm2 :: Floating a => Vector n a -> a
norm2 = sqrt . V.sum . fmap (^2)

-- | ベクトル間のユークリッド距離
distance :: Floating a => Vector n a -> Vector n a -> a
distance xs ys = norm2 (xs ^-^ ys)
```

次に行列と行列演算に関する関数の実装です。

```hs
-- | m × n 行列
type Matrix m n a = Vector m (Vector n a)

-- | リストから行列に変換する
fromList :: (KnownNat m, KnownNat n) => [[a]] -> Matrix m n a
fromList = fromJust . V.fromList . fmap (fromJust . V.fromList)

-- | 全て同じ値を持つ行列を作成する関数
konst :: (KnownNat m, KnownNat n) => a -> Matrix m n a
konst a = V.replicate (V.replicate a)

-- 行列のスカラー倍
(*!!) :: Semiring a => a -> Matrix m n a -> Matrix m n a
a *!! m = mmap (otimes a) m

-- | 各成分ごとの演算
elementWise :: (a -> b -> c) -> Matrix m n a -> Matrix m n b -> Matrix m n c
elementWise op = V.zipWith (V.zipWith op)

-- | 行列の和
(!+!) :: Semiring a => Matrix m n a -> Matrix m n a -> Matrix m n a
(!+!) = elementWise oplus

-- | 行列の差
(!-!) :: Num a => Matrix m n a -> Matrix m n a -> Matrix m n a
(!-!) = elementWise (-)

-- アダマール積
(!!*!!) :: Num a => Matrix m n a -> Matrix m n a -> Matrix m n a
(!!*!!) = elementWise (*)

-- | 内積
dot :: Semiring a => Vector n a -> Vector n a -> a
dot xs ys = V.foldr oplus zero (V.zipWith otimes xs ys)

-- | 行列とベクトルの積
(!*) :: Semiring a => Matrix m n a -> Vector n a -> Vector m a
m !* v = fmap (dot v) m

-- | 転置
transpose :: (KnownNat m, KnownNat n) => Matrix m n a -> Matrix n m a
transpose = V.sequence

-- | 行列積
(!*!) :: (KnownNat p, KnownNat q, KnownNat r, Semiring a)
      => Matrix p q a -> Matrix q r a -> Matrix p r a
a !*! b = flip fmap a \as -> flip fmap b' \bs -> as `dot` bs
  where b' = transpose b

-- 行列の各要素に関数を適用する
mmap :: (a -> b) -> Matrix m n a -> Matrix m n b
mmap f = fmap (fmap f)

-- 行列成分すべての和を取る
msum :: Num a => Matrix m n a -> a
msum = V.sum . fmap V.sum

-- | 隣接行列から最短距離行列を計算する
minDist :: KnownNat m => Matrix m m Double -> Matrix m m Double
minDist = mmap fromTropical . closure . mmap toTropical
  where
    toTropical x = if x == 0.0 then T inf else T x
    fromTropical (T x) = x
```

ずらっと定義を並べましたが一つ一つの実装はとてもシンプルです。

最後の `minDist` は[前述した記事](https://zenn.dev/lotz/articles/9c9ca0708b035b)で実装した `closure` 関数を使って最短距離の行列を作る関数になります。実装はとてもシンプルで、**隣接行列を一度トロピカル代数の世界に飛ばしてから閉包を取り、再び実数の世界に戻してくるだけ**です。

### 最適輸送
グラフのリッチ曲率を定義するためにまずはグラフ上の最適輸送距離を計算できるようにしましょう。`wasserstein`という関数は$xy$成分が$W_1(m_x, m_y)$に対応するような行列を計算します。

```hs
-- | 各ノード間のWasserstein距離を計算する
wasserstein :: KnownNat m
            => Double            -- alpha
            -> Matrix m m Double -- 隣接行列
            -> Matrix m m Double -- 最短距離行列
            -> Matrix m m Double -- Wasserstein距離の行列
wasserstein alpha adj minDist =
  let w = flip V.imap minDist \i ds ->
            flip V.imap ds \j d ->
              if i == j then 0.0 else if i > j then w `V.index` j `V.index` i
                else msum $ mmap (\x -> if isNaN x then 0.0 else x) $ minDist !!*!! shinkhorn (m i) (m j) k
   in w
  where
    eps = 0.2  -- エントロピーによる正則化項の係数
    k   = mmap (\d -> exp (-d / eps)) minDist
    distribute x x' = exp (-minDist `V.index` x `V.index` x')
    m x = let c = V.sum $ V.imap (\x' a -> if a == 0.0 then 0.0 else distribute x x') (adj `V.index` x)
           in flip V.imap (adj `V.index` x) $ \x' a ->
                 if x == x' then alpha
                 else if a == 0.0 then 0.0
                 else (1 - alpha) * distribute x x' / c
```

Wasserstein距離は"距離"なので対角成分は0になり、対称行列になることを利用して必要な計算を半分以下に押さえています。対称行列の性質を利用するところではHaskellの遅延評価を使って、今まさに定義している行列の成分を使うように書くことができ記述を簡潔にできました。最短距離と`shinkhorn`の結果の積は`Infinity * 0.0 = NaN` という計算で`NaN`が出てくることがあるので、`NaN`を`0.0`に変換しています。

`shinknorh`はShinkhornアルゴリズムで最適輸送距離を近似的に計算する関数です。

```hs
-- | Shinkhornアルゴリズム
shinkhorn :: KnownNat m
          => Vector m Double   -- alpha
          -> Vector m Double   -- beta
          -> Matrix m m Double -- K
          -> Matrix m m Double -- P
shinkhorn alpha beta k =
  let delta = 1e-3  -- 収束判定の閾値
      v0 = V.replicate 1.0
      u0 = V.zipWith (/) alpha (k !* v0)
      (u, v) = flip fix (u0, v0) \loop (u, v) ->  -- 収束するまで繰り返す
          let v' = V.zipWith (/) beta (transpose k !* u)
              u' = V.zipWith (/) alpha (k !* v')
           in if (distance v v' + distance u u') < delta then (u', v') else loop (u', v')
   in V.zipWith (*^) u $ transpose (V.zipWith (*^) v (transpose k))
```

実装は上で解説したShinkhornアルゴリズムの式がそのまま対応してる形です。

### グラフのリッチ曲率
グラフのリッチ曲率とリッチフローを実装していきましょう。

```hs
-- | 隣接行列からグラフのリッチ曲率を計算する
ricciCurvature :: KnownNat m => Matrix m m Double -> Matrix m m Double
ricciCurvature adj =
  mmap clean $ (1.0 / (1.0 - alpha)) *!! (konst 1.0 ^-^ elementWise (/) (wasserstein alpha adj md) md)
  where
    md = minDist adj
    alpha = 0.99
    clean x = if x == -inf || isNaN x then 0.0 else x
```

Wasserstein距離行列を最短距離行列で割るとゼロ割りから`Infinity`や`NaN`が発生するので`clean`関数で取り除いています。　

ここまでくればリッチフローの実装はとてもシンプルです。

```hs
-- | 隣接行列からグラフのリッチフローを1ステップ計算する
graphRicciFlow :: KnownNat m => Double -> Matrix m m Double -> Matrix m m Double
graphRicciFlow dt adj = adj - dt *!! ricciCurvature adj !!*!! adj
```

### コミュニティ抽出アルゴリズム
最後にグラフのリッチフローを使ったコミュニティ抽出アルゴリズムを作っていきましょう。

まずはエッジの重みを正規化する関数です。

```hs
-- | 重みの合計がエッジの数に等しくなるように正規化を行う
normalize :: KnownNat m => Matrix m m Double -> Matrix m m Double
normalize m =
  let e = msum (mmap (\x -> if x == 0.0 then 0.0 else 1.0) m)
      c = msum m
   in mmap (\x -> x * e / c) m
```

次にエッジの切断を行う関数です。重みの合計に対して上位5%以内に入る重みを持つエッジを切断します。

```hs
-- | 重みの合計の上位5%を占めるエッジを切断する
surgery :: KnownNat m => Matrix m m Double -> Matrix m m Double
surgery m =
  let total = msum m
      threshold = 0.05 * total
      weights = reverse . sort . concatMap V.toList $ V.toList m
      w_max = flip fix (0, weights) \loop (accum, w:ws) ->
                if accum + w >= threshold then w else loop (accum + w, ws)
   in mmap (\x -> if x >= w_max then 0.0 else x) m
```

これまで実装した関数を組み合わせてコミュニティ抽出のアルゴリズムを実装します。リッチフローは隣接行列の時間発展を計算するので、各ステップでの隣接行列を持つようなリストを計算結果に持つような関数として実装します。コミュニティはリストの最後尾にある隣接行列から連結成分を取り出すことで得られます。またここでの実装はリッチ曲率が収束するまで待つのではなく、予め決められた回数の計算を行うようにしています。

```hs
-- | グラフのリッチフローを利用したコミュニティ抽出
communityDetection :: KnownNat m
                   => Double              -- 微小時間
                   -> Int                 -- イテレーションの回数
                   -> Int                 -- 手術を行うステップ数
                   -> Matrix m m Double   -- 隣接行列
                   -> [Matrix m m Double] -- アルゴリズムによる隣接行列の時間発展
communityDetection eps total s adj = go eps total s adj
  where
    go _   0 _ _ = []
    go eps n s adj =
      let adj'  = graphRicciFlow eps (normalize adj)
          adj'' = if n /= total && n `mod` s == 0 then surgery adj' else adj'
       in adj'' : go eps (n-1) s adj''
```

グラフを使った実験
--------------
最後に作った`communityDetection`を使って実際にグラフのコミュニティ抽出ができているのか見ていくことにしましょう。

隣接行列からグラフを可視化する方法については以下の記事に書いたので良かったら参照してみてください

https://zenn.dev/lotz/articles/9e47c54994876f25d1e5

### コミュニティグラフ
まずは人工的なデータで実験してみましょう。グラフのリッチ曲率の計算例で上げたコミュニティを表すグラフを思い出してください。

![](https://storage.googleapis.com/zenn-user-upload/7m8pgzgv9fje0254vgpnlcwnkjk5 =300x)

このグラフの隣接行列は以下のように書けます。

```hs
adj = fromList [ [0, 1, 1, 1, 1, 0, 0, 0, 0, 0]
               , [1, 0, 1, 1, 1, 0, 0, 0, 0, 0]
               , [1, 1, 0, 1, 1, 0, 0, 0, 0, 0]
               , [1, 1, 1, 0, 1, 0, 0, 0, 0, 0]
               , [1, 1, 1, 1, 0, 1, 0, 0, 0, 0]
               , [0, 0, 0, 0, 1, 0, 1, 1, 1, 1]
               , [0, 0, 0, 0, 0, 1, 0, 1, 1, 1]
               , [0, 0, 0, 0, 0, 1, 1, 0, 1, 1]
               , [0, 0, 0, 0, 0, 1, 1, 1, 0, 1]
               , [0, 0, 0, 0, 0, 1, 1, 1, 1, 0]
               ] :: Matrix 10 10 Double
```

これを使って`communityDetection`による時間発展を可視化したものが以下のGIFアニメーションです。　

![](https://storage.googleapis.com/zenn-user-upload/3df5jllz7z8qvg2idyze7xfl623d =600x)

コミュニティ間をつなぐエッジがどんどん伸びていき、手術によって切断されているのが分かります。真ん中のエッジが切断された後は全ての曲率が同じ値なので安定した配置に収束しています。

### Zachary's karate club
次はより現実的なデータで実験してみましょう。[Zachary's karate club](https://en.wikipedia.org/wiki/Zachary%27s_karate_club) はある大学における空手クラブにおける34人のメンバー間の交友関係をグラフ化したものです。ノードは人を、エッジの存在は交友関係があることを表しています。このグラフのリッチフローによる時間発展をみてみましょう。

![](https://storage.googleapis.com/zenn-user-upload/g05g7n0fsmip9h8bmfhp5jswqrsf =800x)

コミュニティと思える構造がどんどん分離されていく様子が見て取れますね。本当はModularityを使って定量的に評価すべきですが今回は動きを見て満足です。

おわりに
------
ポアンカレ予想に始まり微分幾何・最適輸送、最後にはグラフ理論を少し通ってリッチフローによるコミュニティ抽出のアルゴリズムを眺めてきました。この記事は解説とは言っても必要な事実を並べただけになってしまったので、参考に貼ったリンクからでもそれぞれの理論を深く学ぶキッカケになってもらえれば嬉しいです。どの理論もお話だけではなく緻密に作られた数学の概念がとても綺麗で奥深い世界を作り出しています。そんな抽象的な世界の話をグラフという現実のデータを記述する対象に上手に持ってくることでコミュニティ抽出という役に立つアルゴリズムに使えるというのはとても面白い話ですね。

これは [FOLIO Advent Calendar 2020](https://adventar.org/calendars/5553) の23日目の記事でした。

[^1]: この論文は今年（2020年）の初め頃、まだオフラインでHaskellもくもく会をやっていた頃、一緒にもくもくしていた方に紹介してもらいました
[^2]: 「ある点からの接ベクトルを使った移動」は正確には[指数写像](https://en.wikipedia.org/wiki/Exponential_map_(Riemannian_geometry))によるものを考えます
[^3]: ここでの"距離"は[測地距離](https://ja.wikipedia.org/wiki/%E6%B8%AC%E5%9C%B0%E7%B7%9A)によって測ったものです
[^4]: Ollivier, Yann. "A survey of Ricci curvature for metric spaces and Markov chains." Probabilistic approach to geometry. Mathematical Society of Japan, 2010.
[^6]: Iacovlenco, Olga. "An Introduction to Hamilton’s Ricci Flow."
[^7]: [一意化定理 - Wikipedia](https://ja.wikipedia.org/wiki/%E4%B8%80%E6%84%8F%E5%8C%96%E5%AE%9A%E7%90%86)
[^8]: Ollivier, Yann. "Ricci curvature of Markov chains on metric spaces." Journal of Functional Analysis 256.3 (2009): 810-864.
[^9]: 測度距離空間におけるリッチ曲率の話題としては曲率次元条件を使って下限を評価する話もありますが、それとは別の話です（ですがOllivierのRicci曲率を使ってBonnet-Myersの定理の類似物を示せたりします[^8]）。
[^10]: 実はOllivier-Ricci曲率とは別にCW複体の観点からグラフ上のリッチ曲率の定義を行ったForman-Ricci曲率[^11]というものも提案されていて、こちらの方が組み合わせ的な観点に基づいているため計算自体は高速で行うことができます。両者は相関がありつつもグラフにおける異なる側面を捉えた量であると理解されています[^12]。
[^11]: Sreejith, R. P., et al. "Forman curvature for complex networks." Journal of Statistical Mechanics: Theory and Experiment 2016.6 (2016): 063206.
[^12]: Samal, Areejit, et al. "Comparative analysis of two discretizations of Ricci curvature for complex networks." Scientific reports 8.1 (2018): 1-16.
[^13]: 本来は確率分布ではなく測度という言葉を用いて話すのが適切ですが、ここでは分かりやすさを優先して確率分布という言葉を使って説明していきます。有限集合しか考えないため
[^14]: Lin, Yong, Linyuan Lu, and Shing-Tung Yau. "Ricci curvature of graphs." Tohoku Mathematical Journal, Second Series 63.4 (2011): 605-627.
[^15]: Bourne, David P., et al. "Ollivier--Ricci Idleness Functions of Graphs." SIAM Journal on Discrete Mathematics 32.2 (2018): 1408-1424.
[^16]: Ni, Chien-Chun, et al. "Community detection on networks with ricci flow." Scientific reports 9.1 (2019): 1-12.
[^17]: ここでの正規化の処理はアルゴリズム的な記述になっていますが、微分方程式を修正することで反映することも可能です[^18]。ちなみに元のリーマン多様体のリッチフローでも正規化を考えることができます。
[^18]: Bai, Shuliang, et al. "Discovering hierarchical structures in weighted graphs using Ricci-flow method." arXiv preprint arXiv:2010.01802 (2020).
[^19]: Pythonによる実装は論文の著者らによるものがあります [saibalmars/GraphRicciCurvature](https://github.com/saibalmars/GraphRicciCurvature)
[^20]: データはこちらからダウンロードすることができます。 http://vlado.fmf.uni-lj.si/pub/networks/data/Ucinet/UciData.htm#zachary
[^21]: ここで分類される8つの幾何構造を可視化したサイトもあります <http://www.3-dimensional.space/>