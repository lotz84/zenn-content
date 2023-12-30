---
title: "ベクトルからリストを作る方法 〜次数付きモナドのカン拡張〜"
emoji: "⛓️"
type: "tech"
topics: ["haskell", "圏論", "数学", "math"]
published: true
---
![](https://storage.googleapis.com/zenn-user-upload/yxh1k5rs17yjgs1i4ma8rq3ps7ba)

----

## ベクトルとリスト

**要素を並べたデータ構造**を考える時、

- **ベクトル**は長さが予め（型レベルで）決められたもの
- **リスト**は任意の長さを取れるもの

と区別することがあります。


Haskellの型で表すと、

```hs
data Vector n a = ...

data [] a = [] | a : [a]
```

このように`Vector`は型引数として長さ`n`を取り、リストは型レベルで長さの情報は持っておらず再帰的な型として定義されます。

リストはよく知られているように**モナドの構造**を持っていますがベクトルはどうでしょうか？

実はベクトルにもモナドの構造を与えることができます。

```hs
join :: Vector n (Vector n a) -> Vector n a
```

`join`として"対角線"の要素をとってくる関数を考えれば実際にモナドを作ることができます。

ただこの実装は`Vector n a`を`Finite n -> a`、すなわちnまでのインデックスを与えると要素を返す"関数"と見なしたときのモナド（つまりReaderモナドと同じ挙動）であって、リストモナドの挙動とは少し性質が違います。

リストモナドの挙動に対応するような実装を考えると

```hs
join' :: Vector n (Vector m a) -> Vector (n * m) a
```

このような関数を考えたくなりますが、`n`, `m`という単なるモナドとしては余計な型変数の操作も考慮する必要があるため簡単には行きません。

実はベクトルは**次数付きモナド**[^1]という通常のモナドをちょっと拡張したモナドのインスタンスにすることができます[^2]。

```hs
class GradedMonad (m :: k -> Type -> Type) where
  type Unit m :: k
  type Plus m (i :: k) (j :: k) :: k
  type Inv  m (i :: k) (j :: k) :: Constraint

  pure' :: a -> m (Unit m) a
  bind' :: Inv m i j => m i a -> (a -> m j b) -> m (Plus m i j) b
```

`GradedMonad`のインスタンスとなる型のカインドは`k -> Type -> Type`となっていて、通常のモナドになる型のカインド`Type -> Type`に比べて一つ型引数が多いことが分かります。`GradedMonad`はこの新しい型引数が**型レベルでモノイドの構造を持っている**ことを期待しており、

- `Unit m`が単位元
- `Plus m i j`がモノイドの二項演算

に対応しています。

`pure'`, `bind'`は通常のモナドが持つ`pure`, `>>=`に対応する関数です。

`pure'`は任意の型`a`の値をモノイドの単位元`Unit m`に対応する関手`m (Unit m)`に包まれた型`m (Unit m) a`に写す関数となっています。

また、`bind'`は`Inv m`で与えられる制約を満たす`i`, `j`についてモナドの

```hs
(>>=) :: m a -> (a -> m b) -> m b
```

と同じような挙動をする関数です。この時、カインド`k`の型引数`i`, `j`はモノイドの二項演算`Plus m`によって`Plus m i j`に送られます。

実際にベクトルを`GradedMonad`のインスタンスにしてみましょう。

```hs
import qualified Data.Vector.Sized as V

instance GradedMonad Vector where
  type Unit Vector = 1
  type Plus Vector m n = m * n
  type Inv  Vector m n = ()

  pure' = V.singleton
  bind' = flip V.concatMap
```

`singleton`と`concatMap`は[vector-sized](https://hackage.haskell.org/package/vector-sized)で提供されている関数で以下のような型になっています。

```hs
singleton :: a -> Vector 1 a
concatMap :: (a -> Vector m b) -> Vector n a -> Vector (n * m) b
```

## 次数付きモナド

さて、少し抽象的な話をしましょう。

次数付きモナドはモノイダル圏 $({\mathscr M}, \otimes, I)$ から自己関手のモノイダル圏 $([{\mathscr C}, {\mathscr C}], \circ, {\rm id}_{\mathscr C})$ へのlax monoidal functorです[^3]。

モノイダル圏は既知としましょう（[モノイド圏 - Wikipedia](https://ja.wikipedia.org/wiki/%E3%83%A2%E3%83%8E%E3%82%A4%E3%83%89%E5%9C%8F)）
lax monoidal functor を知るためにまず functor、つまり関手$F$の定義を思い出します。

関手$F$は圏${\mathscr C}$の対象と射を圏${\mathscr D}$の対象と射に対応させ、恒等射と射の合成を保つようなものでした。

今2つのモノイダル圏$({\mathscr C}, \otimes_{\mathscr C}, 1_{\mathscr C})$, $({\mathscr D}, \otimes_{\mathscr D}, 1_{\mathscr D})$を考えると通常の関手であって、更にモノイダル圏の構造を保つような関手を考えることができるでしょう。それは単位対象を単位対象に写し、テンソル積を保つような関手です（あと結合子や単位子とも整合的である必要があります）。

$$
\begin{matrix}
1_{\mathscr D} &=& F(1_{\mathscr C}) \\
F(x)\otimes_{\mathscr D}F(y) &=& F(x \otimes_{\mathscr C} y)
\end{matrix}
$$

このような関手はstrict monoidal functorと呼ばれます。

イコールで結ばれたこれらの条件を単に射が存在するという条件に"緩めた"ものが lax monoidal functor と呼ばれるものです。すなわち以下のような射 $\epsilon$ と自然変換 $\mu$ が存在することを仮定します（あと結合子や単位子とも整合的である必要があります）。

$$
\begin{matrix}
\epsilon&:& 1_{\mathscr D} &\rightarrow& F(1_{\mathscr C}) \\
\mu_{x, y}&:& F(x)\otimes_{\mathscr D}F(y) &\rightarrow& F(x \otimes_{\mathscr C} y)
\end{matrix}
$$

さて、lax monoidal functorが何者であるかが分かったところで、モノイダル圏 $({\mathscr M}, \otimes, I)$ から自己関手のモノイダル圏 $([{\mathscr C}, {\mathscr C}], \circ, {\rm id}_{\mathscr C})$ への lax monoidal functor における$\epsilon$と$\mu$を考えてみましょう。

$$
\begin{matrix}
\epsilon&:& {\rm id}_{\mathscr C} &\rightarrow& F(I) \\
\mu_{i, j}&:& F(i)\circ F(j) &\rightarrow& F(i \otimes j)
\end{matrix}
$$

これらは恒等関手から$F(I)$への自然変換と$F(i)\circ F(j)$（$\circ$ は関手の合成であることに注意）から$F(i \otimes j)$への自然変換なのでHaskellの型として書くと

```hs
pure' :: a -> m (Unit m) a
join' :: m i (m j a) -> m (Plus i j) a
```

のようになり`GradedMonad`の実装と対応していることが分かります。

実は通常のモナドも次数付きモナドとみなすことができます。今、ただ一つの対象からなる離散圏${\mathbb 1}$を自明なモノイダル圏とみなし自己関手のモノイダル圏 $([{\mathscr C}, {\mathscr C}], \circ, {\rm id}_{\mathscr C})$ への lax monoidal functor を考えると、この関手で${\mathbb 1}$における唯一の対象を写した先である自己関手はモナドになります（結合子と単位子との整合性はモナド則に対応します）。

ところで、モノイダル圏${\mathscr M}$から自己関手の圏$[{\mathscr C}, {\mathscr C}]$への関手 $H$ として次数付きモナドがあって、モノイダル圏${\mathbb 1}$から自己関手の圏$[{\mathscr C}, {\mathscr C}]$への関手 $L$ としてモナドが考えられるなら、次数付きモナド$H$の左カン拡張としてモナド $L$ が得られるのか考えたくなりますよね？

![](https://storage.googleapis.com/zenn-user-upload/mucbphykq2aeykaj3kjpb2ipy0kq =400x)

特に今はベクトルという次数付きモナドとリストというモナドを考えているので、もし**ベクトルのカン拡張としてリストが得られる**としたらとても面白そうですよね。

## ベクトルからリストへ

これを実際にHaskellで実装して確かめてみましょう。Haskellではコエンドを利用して以下のように左カン拡張を定義することができます。

```hs
data Lan g h a where
  Lan :: (g b -> a) -> h b -> Lan g h a
```

`Lan g h a`は`h`の`g`に沿った左カン拡張になっています。`b`が存在型として隠蔽されているのがポイントです。この実装は[kan-extentions](http://hackage.haskell.org/package/kan-extensions-5.2.1/docs/Data-Functor-Kan-Lan.html)というライブラリで提供されています。

これを使ってベクトルのカン拡張としてリストを定義してみましょう。

```hs
newtype Flip f a b = Flip { runFlip :: f b a }

newtype List a = List { getList :: Lan (Const ()) (Flip Vector a) () }
```

`Flip` は `Vector` の型引数の順番を入れ替えるための型です。`List a` は `Flip Vector a` の `Const ()` に沿った左カン拡張として定義されています。

このように定義した`List a`が実際に`[a]`と同型であることを確認してみましょう。

```hs
-- | ベクトルをリストとみなすための補助関数
fromVector :: Vector n a -> List a
fromVector v = List (Lan (const ()) (Flip v))

toList :: List a -> [a]
toList (List (Lan _ (Flip v))) = V.toList v

fromList :: [a] -> List a
fromList xs = V.withSizedList xs (\v -> fromVector v)
```

これらの関数によって`List a`は情報を落とすこと無く通常のリストと相互に変換できることが分かります。

実は`List a`の型をよく見て変形していくと

```hs
Lan (Const ()) (Flip Vector a) ()

<=> {- Lan と Flip の定義より -}

exists n. (Const () n -> (), Vector n a)

<=> {- Const () n <=> (), () -> () <=> (), ((), A) <=> A より -}

exists n. Vector n a
```

となり、単に**リストとはなんらかの長さを持ったベクトルである**ということを表しています。ただし `exists n.`は`n`が存在型であることを表しています（残念ながらHaskellの文法には直接は存在しない記号です）。

さて、ここまでは**データ構造として**ベクトルのカン拡張がリストになる事を見てきましたが、ここからはこの構成によってベクトルの次数付きモナドの性質から**リストモナドの性質をどこまで復元できるか**見ていきましょう。

まず、FunctorとApplicativeのインスタンスは以下のように実装することができます。

```hs
instance Functor List where
  fmap f (List (Lan _ (Flip v))) = fromVector $ fmap f v

instance Applicative List where
  pure a = fromVector $ pure' a
  (List (Lan _ (Flip vf))) <*> (List (Lan _ (Flip va))) =
    fromVector $ vf `bind'` \f -> fmap f va
```

面白いことに`GradedMonad`の実装から`Functor`と`Applicative`の実装は自動的に手に入ってしまうのです（つまりこれはベクトルに限った話ではありません）。実行してみると期待通りの実装になっていることが分かります。

```hs
> toList $ fmap (+1) $ fromList [1,2,3]
[2,3,4]

> toList $ (fromList [(+1), (*2)]) <*> (fromList [1,2])
[2,3,2,4]
```

いよいよモナドを実装しましょう。しかし少し試すと分かりますが、これは簡単には行きません。一般に次数付きモナドのカン拡張が存在したとしてもそれがまた次数付きモナドになっているとは限らないのです。["A Criterion for Kan Extensions of Lax Monoidal Functors"](https://arxiv.org/abs/1809.10481)には（強モノイダル関手に沿った）次数付きモナドが再び次数付きモナドになるための十分条件（Theorem 2.1）が書かれています。この条件を今回のベクトルとリストの例に翻訳すると大まかには、

`List (List a)` と `Lan (Const ()) (\(n, m) -> Vector n (Vector m a)) ()`

の間に同型対応が存在することを要求しています（一部Haskellに翻訳しきれていない部分がありますが察して下さい）。この条件は`List (List a)`の値があるn×mの二次元配列で表せることを要求しており非常に強く、一般的には成り立ちません。ここではこの条件が成り立てば成り立つという意味で、より弱い関数の存在を仮定しましょう。

```hs
joinL :: List (List a) -> List a
joinL = fromList . concatMap toList . toList
```

ほぼ欲しかったものの存在を仮定しちゃいましたね😅 ただこれは`List (List a)`の値がある長さのベクトルで表せるという仮定なので先程の条件よりは弱いものになっています。これを使えば、

```hs
instance Monad List where
  m >>= k = joinL $ fmap k m
```

とめでたくモナドの実装を行うことができました👏

## あとがき
この話は[graded monad in nLab: 3. Uses](https://ncatlab.org/nlab/show/graded+monad#uses)を読んで知りました。Haskellでどこまで再現できるか分かりませんでしたが思ったよりうまく行ってホッとしています。無事、次数付きモナドとみなしたベクトルの左カン拡張としてリストモナドの実装を手に入れることができました。最後モナドの実装だけ強い関数の存在を仮定してしまいましたが、もしかすると緩い条件で十分かもしれません（もし知っていたらこっそり教えて下さい）。今回と同様の構成で他にもよく知られたHaskellのデータ型を`GradedMonad`の`Lan`として実装できるものがあるかもしれません。データ構造の間にある深い関係が分かってくるとワクワクしますよね。

最後にこの話を書くために深夜"2時"まで議論に付き合ってくれたryota-ka氏に感謝🙏

----

＼読んでいただきありがとうございました！／
この記事が面白かったら いいね♡ をいただけると嬉しいです☺️
バッジを贈っていただければ次の記事を書くため励みになります🙌

[^1]: "次数付き"という名前は次数付き環に由来しています。
Fujii, Soichiro, Shin-ya Katsumata, and Paul-André Melliès. "Towards a formal theory of graded monads." International Conference on Foundations of Software Science and Computation Structures. Springer, Berlin, Heidelberg, 2016.
[^2]: ここでの実装は Dominic Orchard, Tomas Petricek, Embedding effect systems in Haskell, ([pdf](https://www.doc.ic.ac.uk/~dorchard/publ/haskell14-effects.pdf)) を参考にしています。
[^3]: 何か問題でも？