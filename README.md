# Project 12-Covid19

執行程式project12.m後，
會呈現出一張感染人數、累積康復人數隨時間變化的圖表，還有一張R0值隨時間變化的的折線圖，在command window也會以文字呈現每天的R0值；
由於gamma值有可能為零或趨近於零，R0會因此不存在或極大，使0~10之間的數值變化難以直接用圖表觀察;因此不存在的話則在圖表上數值-20的位置標記橘色的x ，在R0超過20的情況下皆以20呈現且同樣在圖表上數值-20的位置標記橘色的 x
（如果R0值存在，橘色x 則標記在數值0的位置）


# SIR Model
S.I.R. 分別是 **Susceptible (易感染)、Infected (已感染)、Removal (移除感染)**
的縮寫， S.I.R. 模型的建構需有幾個基本假設，
* **第一**，不考慮個體的出生、自然死亡與流動的情況，個體總數是固定數**N**；
* **第二**，每一個體與其他個體接觸的機會均等。

假設某區域的易感染人數、已感染人數、移除感染人數均是時間 t 的函數，
分別記作 S(t), I(t), R(t), 得 S(t) + I(t) + R(t) = N .
**β**(Beta) 值的意義是一個已感染者每天有效傳播病毒的人數占總人口的比例，而 **γ**(Gamma) 為每天移除感染的人數占已感染人數的比例

## **R0** :Basic reproduction number
為一個病人在疫情期間直接感染他人個數的平均值: **R0=N*Beta./Gamma**
* (the average number of persons directly infected by an infectious case during the entire infectious period 
when s/he enters a totally susceptible population)
---
* dSdt=−(β⋅S⋅I)/N
* dIdt=[(β⋅S⋅I)/N]−ν⋅I
* dRdt=ν⋅I

既然已知 S+I+R=N，將S的值替換為 (N-R-I) 
我們可以將式子改寫成:

* dIdt = [βI(N-I-R)/N]-vI  
* dRdt= vI

若已知一個國家每日的感染人數(I)及康復人數(R)，設將時間間隔dt設為1天，dI,dR分別是感染人數和康復人數的變化量，tt為目前日期

* dI=(beta(tt+1).* I(tt).*(N-I(tt)-R(tt)) - gamma(tt+1).* I(tt)).* dt
* dR=(gamma(tt+1).* I(tt)).*dt

則R(tt)+dR為明天的累積康復人數，I(tt)+dI為明天的目前感染總人數
，其中beta(tt+1),gamma(tt+1)為未知數

藉由**fsolve**

* fsolve(@(x) R(tt+1)-(x.* I(tt)).* dt-R(tt),1) 得到的 **x**值為 gamma(tt+1)
* fsolve(@(y) I(tt+1)-(y.* I(tt).* (N-I(tt)-R(tt)) - gamma(tt+1).* I(tt)).* dt-I(tt),1) 得到的 **y**值為beta(tt+1)

既然有了每天的beta和gamma值，即可求得每日的**R0**值



---




* 台灣2020/01/24-2020/12/14累積感染人數＆累積康復人數 資料來源：[Link](https://zh.m.wikipedia.org/wiki/Template:2019%E5%86%A0%E7%8B%80%E7%97%85%E6%AF%92%E7%97%85%E7%97%85%E4%BE%8B%E6%95%B8/%E8%87%BA%E7%81%A3%E5%9C%96%E8%A1%A8)
