## 3.4.3 平均・精度が未知である場合

ガウス・ウィシャート分布を用いれば事後分布が同じ形式になり、便利である。

 $$
 \begin{eqnarray}
 p(\boldsymbol{\mu,\Lambda}) &=& NW(\boldsymbol{\mu,\Lambda}|\boldsymbol{m},\beta,\nu,\boldsymbol{W}) \\
 &=& \mathcal{N}(\boldsymbol{\mu}|\boldsymbol{m},(\beta\boldsymbol{\Lambda})^{-1})\mathcal{W}(\boldsymbol{\Lambda}|\nu,\boldsymbol{W})\tag{3.125}
 \end{eqnarray}
 $$

 1次元データの場合の式 $(3.81)$ と同様に $(3.125)$ も右から左に見ると、 $\boldsymbol{\Lambda}$ が分かったら $\boldsymbol{\mu}$ が出る と見ることが可能で、左から右に見た場合は、最初は $\boldsymbol{\Lambda}$ は既知として良い、という気持ちで式を眺めると。

 $$p(\mu, \lambda) = \mathcal{N}(\mu|m,{(\beta\lambda)}^{-1})\mathrm{Gam}(\lambda|a,b) \tag{3.81}$$

 最初のfix parameterは一次元の場合と同様に $\boldsymbol{m},\beta,\nu,\boldsymbol{W}$ であり、これらのパラメータをデータを（あるいはデータの集合）を見るたびに更新していく。

 以前と全く同様に、$(3.128)$ 式では $\boldsymbol{\Lambda}$ は一旦既知として、データを観測したあとに $\boldsymbol{m}\to\boldsymbol{\hat{m}},\beta\to\hat{\beta}$ と $\boldsymbol{\mu}$ のパラメータを更新する。1次元の場合と全く同じ流れで、 $(3.97),(3.102),(3.103)$ を以下のように書き換える。
$$
\begin{eqnarray}
p(\boldsymbol{\mu}) &=& \mathcal{N}(\boldsymbol{\mu}|\boldsymbol{m},(\boldsymbol{\Lambda_{\boldsymbol{\mu}}})^{-1}) \tag{3.97} \\
p(\boldsymbol{\mu}) &=& \mathcal{N}(\boldsymbol{\mu}|\boldsymbol{m},(\beta\boldsymbol{\Lambda})^{-1}) \tag{3.97'} \\ \\
\boldsymbol{\hat{\Lambda}_{\boldsymbol{\mu}}} &=& N\boldsymbol{\Lambda} + \boldsymbol{\Lambda_{\boldsymbol{\mu}}} \tag{3.102} \\
\hat{\beta}\boldsymbol{\Lambda} &=& N\boldsymbol{\Lambda} + \beta\boldsymbol{\Lambda} \tag{3.102'} \\ \\
\hat{\boldsymbol{m}} &=& \boldsymbol{\hat{\Lambda}_{\boldsymbol{\mu}}}^{-1} (\boldsymbol{\Lambda}\sum_{i=0}^N \boldsymbol{x}_n + \boldsymbol{\Lambda_{\boldsymbol{\mu}}}\boldsymbol{m}) \tag{3.102} \\

\hat{\boldsymbol{m}} &=& {(\hat{\beta}\boldsymbol{\Lambda})}^{-1} (\boldsymbol{\Lambda}\sum_{i=0}^N \boldsymbol{x}_n + \beta\boldsymbol{\Lambda}\boldsymbol{m}) \tag{3.102}

\end{eqnarray}
$$  


一次元の場合と同じ流れで、 $\boldsymbol{\Lambda}$ の事後分布 $p(\boldsymbol{\lambda}|\boldsymbol{X})$ を求める。
一次元の場合、p.95では $p(\lambda|\boldsymbol{X})$ を求めるために
$$p(\boldsymbol{X},\mu,\lambda) = p(\mu|\lambda,\boldsymbol{X})p(\lambda|\boldsymbol{X})p(\boldsymbol{X}) \tag{3.84} $$

を用いて
$$
\begin{eqnarray}
p(\lambda|\boldsymbol{X}) &=& \frac{p(\boldsymbol{X},\mu,\lambda)}{p(\mu|\lambda,\boldsymbol{X})p(\boldsymbol{X})} \\ \\ &\propto& \frac{p(\boldsymbol{X},\mu,\lambda)}{p(\mu|\lambda,\boldsymbol{X})}
\\ \\
&=& \frac{p(\boldsymbol{X}|\mu,\lambda)p(\mu,\lambda)}{p(\mu|\lambda,\boldsymbol{X})}
\tag{3.85'}\end{eqnarray} $$

と表せる。 $\mu\to\boldsymbol{\mu},\lambda\to\boldsymbol{\Lambda}$ として
$$ p(\boldsymbol{\Lambda}|\boldsymbol{X}) \propto \frac{p(\boldsymbol{X}|\boldsymbol{\mu,\Lambda})p(\boldsymbol{\mu,\Lambda})}{p(\boldsymbol{\mu}|\boldsymbol{\Lambda},\boldsymbol{X})}$$

の対数を取り（比例）定数を無視すれば、

$$
\begin{eqnarray}
p(\boldsymbol{\Lambda}|\boldsymbol{X}) &=& \sum \mathrm{ln}\mathcal{N}({\boldsymbol{x}}_{n}|\boldsymbol{\mu,\Lambda^{-1}}) \tag{1}\\ \\
&+&  \mathrm{ln}\mathcal{N}({\boldsymbol{\mu}}|\boldsymbol{m,{(\beta\Lambda)}^{-1}}) \tag{2} \\ \\
&+& \mathrm{ln}\mathcal{W}(\boldsymbol{\Lambda}|\nu,\boldsymbol{W}) \tag{3} \\ \\
&-& \mathrm{ln}\mathcal{N}({\boldsymbol{{\mu}}}|\boldsymbol{\hat{m},{(\hat{\beta}\Lambda)}^{-1}}) \tag{4}
\end{eqnarray}$$

$(1),(2),(3),(4)$ の最終結果だけ書くと、
$$\begin{eqnarray}
(1) &=& \frac{N}{2}\mathrm{ln}|\boldsymbol{\Lambda}| -\frac{1}{2}\sum \mathrm{Tr}(({{\boldsymbol{x}}_n}{{\boldsymbol{x}}}^{\top}_n -
  {{\boldsymbol{x}}_n}{{\boldsymbol{\mu}}}^{\top} -
  {{\boldsymbol{\mu}}}{{\boldsymbol{x}}}^{\top}_n +
  {{\boldsymbol{\mu}}}{{\boldsymbol{\mu}}}^{\top} )\boldsymbol{\Lambda}) \\ \\

(2) &=& \frac{1}{2}\mathrm{ln}|\beta\boldsymbol{\Lambda}|-\frac{\beta}{2}\mathrm{Tr}(({{\boldsymbol{\mu}}}{{\boldsymbol{\mu}}}^{\top} -
  {{\boldsymbol{\mu}}}{{\boldsymbol{m}}}^{\top} -
  {{\boldsymbol{m}}}{{\boldsymbol{\mu}}}^{\top} +
  {{\boldsymbol{m}}}{{\boldsymbol{m}}}^{\top} )\boldsymbol{\Lambda}) \\ \\

(3) &=& \frac{\nu-D-1}{2}\mathrm{ln}|\boldsymbol{\Lambda}| - \frac{1}{2}\mathrm{Tr}(\boldsymbol{W}^{-1}\boldsymbol{\Lambda}) \\ \\
(4) &=& -\frac{1}{2}\mathrm{ln}|\hat{\beta}\boldsymbol{\Lambda}| + \frac{\hat{\beta}}{2}
\mathrm{Tr}(({{\boldsymbol{\mu}}}{{\boldsymbol{\mu}}}^{\top} -
  {{\boldsymbol{\mu}}}{{\boldsymbol{\hat{m}}}}^{\top} -
  {{\boldsymbol{\hat{m}}}}{{\boldsymbol{\mu}}}^{\top} +
  {{\boldsymbol{\hat{m}}}}{{\boldsymbol{\hat{m}}}}^{\top} )\boldsymbol{\Lambda})
\end{eqnarray}$$
となる。 対称行列は対角化が可能であることにより、トレースを用いた。$(3.129)$ を変形し

$$\begin{eqnarray}
\hat{\beta} &=& N + \beta \\ \\
\hat{\boldsymbol{m}} \hat{\beta} &=& \sum {\boldsymbol{x}}_n + \beta\boldsymbol{m}
\end{eqnarray}
$$

これを用いて $(1)\sim(4)$ を整理する。途中式はpdf参照。通常のウィシャート分布のlogの形である $(2.87)$ と比較して、
$$\begin{eqnarray}
\mathrm{ln}\mathcal{W}(\boldsymbol{\Lambda}|\nu,\boldsymbol{W})
&=& \frac{\nu-D-1}{2}\mathrm{ln}|\boldsymbol{\Lambda}| - \frac{1}{2}\mathrm{Tr}(\boldsymbol{W}^{-1}\boldsymbol{\Lambda}) +\mathrm{const.} \tag{2.87} \\\\
\mathrm{ln}p(\boldsymbol{\Lambda}|\boldsymbol{X}) &=&\frac{N+\nu-D-1}{2}\mathrm{ln}|\boldsymbol{\Lambda}| \\ \\
&-& \frac{1}{2} \mathrm{Tr}((\sum{{\boldsymbol{x}}_n}{{\boldsymbol{x}}}^{\top}_n +\beta
  {{\boldsymbol{m}}}{{\boldsymbol{m}}}^{\top} -\hat{\beta}
  {{\boldsymbol{\hat{m}}}}{{\boldsymbol{\hat{m}}}}^{\top} +
  \boldsymbol{W}^{-1})\boldsymbol{\Lambda}) \tag{3.131} \\ \\ &+& \mathrm{const.}
\end{eqnarray}$$

超パラメータは $(3.133)$ のようにすれば良いと分かる。
$$ \begin{eqnarray}
\boldsymbol{\hat{W}}^{-1} &=& \sum{{\boldsymbol{x}}_n}{{\boldsymbol{x}}}^{\top}_n +\beta
  {{\boldsymbol{m}}}{{\boldsymbol{m}}}^{\top} -\hat{\beta}
  {{\boldsymbol{\hat{m}}}}{{\boldsymbol{\hat{m}}}}^{\top} +
  \boldsymbol{W}^{-1}
  \\ \tag{3.133}\\
  \hat{\nu} &=& N + \nu
\end{eqnarray}$$

予測分布の更新時にはベイズ公式及び巻末のウッドベリーの公式を用いて、ステューデントのt分布と比較して更新パラメータを割り当てれば良い。
