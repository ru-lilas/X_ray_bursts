# Todo
## understand synchrotron radiation --> Rybicki & Lightman
   1. [x] 6.1 
   2. [x] 6.2
   3. [x] 6.3
   4. [ ] 6.8 SSA
## low mass X-ray binary
1. [ ] accretion rateに対するX-ray bursterのpopulationを調べる
   - NinjaSat論文の結果が他のsystemでも成立するのなら, X-ray bursterのaccretion rateはan order - two orderで変わることになる。
    ~~LMXB accretion rate は0.1L_Edd - 0.01L_Eddに多い --> X-ray burstを起こすのはどれか~~
     - [ ] そのためにはまずX-ray bursterたちのluminosityを集める必要がある
     - 'annual review of astronomy and astrophysics'からX-ray bursterのカタログをつくると良さそう？
     - ちょっと古いが[https://www.annualreviews.org/content/journals/10.1146/annurev.astro.44.051905.092519]に情報ありそう
2. [ ] clocked bursterのpropertyをまとめた表をつくる(clockedは6天体しか同定されていない)
   1. [ ] まずは6天体のお名前を知る
   2. [ ] recurrence time
   3. [ ] distance
   4. [ ] peak luminosity --> accretion rateにつながる
   - luminosityに何か規則性はあるのか？
3. [ ] not-clocked X-ray burstの recurrence timeを調べる
   - recurrence timeの感覚を身につける（どれくらいのoederなのか？）
4. [ ] LMXB の典型的な公転周期（そもそも典型的なものがあるかどうか）を調べる
5. [ ] lightcurveやburstのタイムスケールがどんな物理で決まっているのか？

# ejecta
ある半径rにある球対称密度n_eの速度v_shock で進む衝撃波によって加速された電子からのシンクロトロン放射だと仮定する
1. シンクロトロン放射の放射スペクトルを求めよ（モデルを適用）
2. XRBが発生してすぐに0.3cでものが飛んでいる. 3min後にどれくらい進んでいるかというとr_{shock} = v_{wind} \times t = 1.8e12 cm (v/0.3c)(t/3min)
3. 飛んでいったものの質量m_ejが知りたい. 電波が落ちるまでの時間tau_{decay}=10 min でnon-relの衝撃波を減速させるには, どれくらいの質量とぶつかればよいか？ --> だいたい同じくらいの質量とぶつかればよい. 4pi (r_sh)^2 n_{amb}v_{sh}tau = M_{ej} --> M_ejが決まったとすると n_{amb}が決まる --> forward shockでどれくらいn_{amb}の電子が加速され, シンクロトロン放射するか？
4. どれくらいのエネルギーの電子が, どれくらいの時間で加速されるのか？ L_{shock} <= 4pi r_{sh}^2 (n_{namb}m v_{sh}^2 / 2) v_{sh}, releaseされるエネルギーのうち一部は磁場に ~ epsilon_B B^2 / 4pi -- この式からB_{sh}でてくる
5. 電子の分布はpower lawに従うとしよう dN/d gamma ~ C gamma^{-p} (gamma_{min} < gamma < gamma_{max}) L_sh * tau ~ epsilon_e int{dN/d gamma}gamma m c^2 d gamma <-- このしきからC出てくる
6. M_{ej}, p, epsilon_{B}, epsilon_{e}はparameter
   - M_{ej}はrecurrence timeの間に降着してきたmassのeta倍が飛んでいくとしてみる
7. 無衝突衝撃波を仮定したときに得られる電子のエネルギー(Lorentz factor)がわかれば, シンクロトロン放射のspectrumも得られる(ほんとに?)

1. X-ray burstからどのようにしてshockが形成されるか
   - shockはどれくらいの速さ？
2. shockが外側の物質を加熱する
   - 外側ってどこ？
   - envelope?
3. X-ray burstが起きてから電波が立ち上がるまでの時間は~3min
   - shockによる加熱は3min以内に完了しなければならない
4. 加熱されたelectronが相対論的に運動-->シンクロトロン放射
   - shock自身が放射するわけではない？
5. 電波はX-ray burstから~10minで落ちる
   - シンクロトロン放射している電子は~10min程度で冷却されなければならない
