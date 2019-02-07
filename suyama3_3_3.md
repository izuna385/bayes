# 3.3.3　平均・精度が未知の場合


$$ p(x|\mu,\lambda) = \mathcal{N}(x|\mu,{\lambda}^{-1}) \tag{3.80}$$

左辺について、変数の明示の規則は $(3.47),(3.64)$ などと同様に、未知数の場合は明示する。  

　左辺→右辺の気持ちは、 $\mu,\lambda$ が分かっていない場合にデータ $x = x$ が降ってくる確率をガウス分布に問うたなら、ガウス分布の式 $(2,64)$  に $x,\mu,{\lambda}$ を代入すれば確率を吐き出してくれる、と解釈できる。  
(ただし、3.3.3では精度パラメータ ${\lambda}^{-1} \equiv {\sigma}^{2}$ とする。)

 ![](./img/i2.png)  


$m,\beta,a,b$ を固定パラメータとして、ある $\mu = \mu, \lambda = \lambda$ が（同時に）発生する確率  $p(\mu,\lambda)$ 、つまり $\mu, \lambda$ の事前分布として $(3.81)$ を用いると、同じ形式の事後分布を得られることが分かっている。  

 $p(\mu, \lambda) = \mathcal{N}(\mu|m,{(\beta\lambda)}^{-1})\mathrm{Gam}(\lambda|a,b)\tag{3.81}$

 　この章での目標は、新しい観測データ$x$を観測したとき(=事後)、その分布すなわち事後分布を決めるような超パラメータ $\hat{m},\hat{\beta},\hat{a},\hat{b}$ を求めることである。

$(3.81)$ では $\mu,\lambda$ は独立ではないが（図3.5の3番目）、それぞれの以前導出した結果を流用することを考える。   
 
  ここで、$(3.81)$ 式を条件付き確率を見るように眺めると、以降少しだけ気持ちが楽になる。つまり、 $\mathrm{Gam}(\lambda|a,b)$ によって $\lambda$
は既知だと分かったという条件下での $\mathcal{N}(\mu|m,{(\beta\lambda)}^{-1})$  と $(3.81)$ を見ると、 $\mu$ を求めるときに少しだけ楽になる。
## 平均 $\mu$ の更新
初期固定値 $m,{\lambda}_{\mu}$ を適当に定めた場合、平均が未知であるときにガウス分布から $x$ が返ってくる確率は  
$$p(x|\mu) = \mathcal{N}(x|\mu,{\lambda}^{-1})\tag{3.47}$$

 ただし $\mu$ の事前分布として以下のように定める。  

   $$ p(\mu)=\mathcal{N}(\mu|m,{{\lambda}_{\mu}}^{-1})\tag{3.48}$$

  上のように定めた場合に、 $\mathbf{X}={\{x_1,...,x_N\}}$ を観測したときの超パラメータ $\hat{\lambda}_{\mu},\hat{m}$ は以下のように更新される。  
  
  $$\hat{\lambda}_{\mu} = N\lambda + \lambda_\mu\tag{3.53}$$  
  $$\hat{m} = \frac{\lambda\sum_{i = 1}^Nx_n+\lambda_\mu}{\hat{\lambda_\mu}}\tag{3.54}$$

  $(3.81)$ の $\mathcal{N}(\mu|m,{(\beta\lambda)}^{-1})$ と $(3.48)$ の $p(\mu)=\mathcal{N}(\mu|m,{{\lambda}_{\mu}}^{-1})$ とを比べて、 $\lambda_\mu = \beta\lambda$ と置き換えたいが $(3.51)$ と比較して $(3.82)$ の意味するところを考える。

   $$p(\mu|\mathbf{X}) = \mathcal{N}(\mu|\hat{m},\hat{{{\lambda}_{\mu}}}^{-1})\tag{3.51} $$

  $(3.51)$ 導出時は $\lambda$ は既知であったが、今回も気持ちは同じである。つまり、$\lambda$ を更新する超パラメータの考慮等はすべて　$(3.81)$ 式 $Gam$ 以降に押し付けるので、今は $\lambda$ は既知として良い。$\lambda$ は一旦既知でfixであるとして、超パラメータの更新を $\hat{\beta}$ に押し付けることを考えると、$(3.53)$ は次のように書き換えられる。右辺の $\lambda_\mu$ についてはそのまま $\beta\lambda$ として、真ん中の $N\lambda$ 内の $\lambda$ はそのまま（今は一旦既知としている）$\mathbf{X}$ の精度と見て良い。（新しいデータが来たときの $\lambda$ の更新は、 $Gam$ ですべて考えるので）左辺について、 $\lambda$ は一旦既知として、更新分を $\hat{\beta}$ に押し付ける。

   $$\hat{\beta}\lambda = N\lambda + \beta\lambda \tag{3.53'}$$

  $(3.81)$ 上では,
  $$ \mathcal{N}(\mu|m,{(\beta\lambda)}^{-1}) \to \mathcal{N}(\mu|m,{(\hat{\beta}\lambda)}^{-1})$$

   を行っていることに対応する。
## 精度 $\lambda$ の更新
  式 $(3.85)$ を書き下すと以下のようになる。

  $$\begin{eqnarray}p(\lambda|\mathbf{X}) &=& \frac{{\{\prod_{n=1}^{N}\mathcal{N}(x_n|\mu,\lambda^{-1})\}}\mathcal{N}(\mu|m,{(\beta\lambda)}^{-1})\mathrm{Gam}(\lambda|a,b)}{\mathcal{N}(\mu|\hat{m},{(\hat{\beta}\lambda)}^{-1})} \\ &=& \frac{\left[\prod_{n=1}^{N}[\sqrt{\frac{\lambda}{2\pi}}exp(-\frac{\lambda}{2}{(x_i-\mu)}^2)]\right]\left[\sqrt{\frac{\beta\lambda}{2\pi}}exp(-\frac{\beta\lambda}{2}{(\mu-m)}^2)\right]\left[C_G(a,b)\lambda^{a-1}e^{-b\lambda}\right]}{\sqrt{\frac{\hat{\beta}\lambda}{2\pi}}exp(-\frac{\hat{\beta}\lambda}{2}(\mu-\hat{m})^2)}\end{eqnarray}$$

この式に従い, $(3.83)$ 式を適宜代入して計算すると、$\lambda$ の更新を担う超パラメータの更新は教科書 $(3.88)$ となる。（途中式はノート参照）

　結局、ガンマ・ガウス分布を用いることで、(いつもと同じような)超パラメータがデータのみによって更新される形に帰着される。

## $p(x^{ * })$ の導出
 $(3.90)$  の第一項目に $(3.80)$ を、第二項目に、$N=1$ として $(3.83),(3.88)$ によって超パラメータを更新したガウス・ガンマ分布 $(3.91)$ を用いれば良い。
　何も観測していない状態で勝手に超パラメータを決め打ちした場合は、スチューデントのt分布を用いることになるが、実際はデータを $N$ 個観測したあとの予測分布を用いることが多い。
