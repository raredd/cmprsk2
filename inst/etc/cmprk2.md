---
date: "05 Feb 2019"
output:
  html_document:
    keep_md: yes
---

cmprsk2
====

Extensions for the [cmprsk](https://cran.r-project.org/web/packages/cmprsk/index.html) package.

See [intro vignette](vignettes/cmprsk2.Rmd)

To install:

```r
# install.packages('devtools')
devtools::install_github('raredd/cmprsk2', build_vignettes = TRUE)
```



### crr formula method


```r
## model deaths with ltx and withdraw as competing events
cr1 <- crr2(Surv(futime, event(censored) == death) ~ age + sex + abo,
            data = transplant)

## include an all-cause death overall survival model
cr2 <- crr2(Surv(futime, event(censored) == death) ~ age + sex + abo,
            data = transplant,
            cox = Surv(futime, event == 'death') ~ age + sex + abo)
```

### crr2 summary methods


```r
summary(cr1)
```

```
## $`CRR: death`
##         HR  L95  U95    p
## age   1.02 1.00 1.04 0.09
## sexf  0.65 0.40 1.07 0.09
## aboB  1.48 0.70 3.12 0.31
## aboAB 1.14 0.34 3.87 0.83
## aboO  1.48 0.85 2.56 0.17
## 
## $`CRR: ltx`
##         HR  L95  U95    p
## age   0.99 0.99 1.00 0.15
## sexf  1.09 0.93 1.29 0.28
## aboB  0.70 0.54 0.91 0.01
## aboAB 0.91 0.59 1.39 0.65
## aboO  0.59 0.49 0.70 0.00
## 
## $`CRR: withdraw`
##         HR  L95   U95    p
## age   0.98 0.95  1.02 0.28
## sexf  1.42 0.75  2.68 0.28
## aboB  2.35 0.82  6.74 0.11
## aboAB 3.04 0.81 11.43 0.10
## aboO  2.37 1.04  5.41 0.04
```

```r
library('htmlTable')
summary(
  cr2,
  html = TRUE, n = TRUE, ref = TRUE,
  htmlArgs = list(
    caption = 'CRR models.',
    rgroup = c('Age', 'Sex', 'Blood type'),
    rnames = c('+1 year change', 'Female', 'B', 'AB', 'O'),
    css.cell = 'white-space: nowrap; padding: 0px 5px 0px;'
  )
)
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr><td colspan='18' style='text-align: left;'>
CRR models.</td></tr>
<tr>
<th style='border-top: 2px solid grey;'></th>
<th colspan='1' style='font-weight: 900; border-top: 2px solid grey; text-align: center;'></th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
<th colspan='3' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Cox PH</th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
<th colspan='3' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>CRR: death</th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
<th colspan='3' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>CRR: ltx</th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
<th colspan='3' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>CRR: withdraw</th>
</tr>
<tr>
<th style='border-bottom: 1px solid grey;'> </th>
<th style='border-bottom: 1px solid grey; text-align: center;'>Total<br /><font size=1>n = 797 (%)</font></th>
<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>Events<br /><font size=1>n = 66 (8)</font></th>
<th style='border-bottom: 1px solid grey; text-align: center;'>HR (95% CI)</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>p</th>
<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>Events<br /><font size=1>n = 66 (8)</font></th>
<th style='border-bottom: 1px solid grey; text-align: center;'>HR (95% CI)</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>p</th>
<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>Events<br /><font size=1>n = 618 (78)</font></th>
<th style='border-bottom: 1px solid grey; text-align: center;'>HR (95% CI)</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>p</th>
<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>Events<br /><font size=1>n = 37 (5)</font></th>
<th style='border-bottom: 1px solid grey; text-align: center;'>HR (95% CI)</th>
<th style='border-bottom: 1px solid grey; text-align: center;'>p</th>
</tr>
</thead>
<tbody> 
<tr><td colspan='18' style='font-weight: 900;'>Age</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;+1 year change</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>797 (100)</td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>66 (100)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.02 (0.99, 1.04)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.16</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>66 (100)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.02 (1.00, 1.04)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#990000">0.090</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>618 (100)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>0.99 (0.99, 1.00)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.15</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>37 (100)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>0.98 (0.95, 1.02)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.28</font></td>
</tr> 
<tr><td colspan='18' style='font-weight: 900;'>Sex</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Female</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>438 (55)</td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>42 (64)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>42 (64)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>332 (54)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>17 (46)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;B</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>359 (45)</td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>24 (36)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>0.69 (0.42, 1.14)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.15</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>24 (36)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>0.65 (0.40, 1.07)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#990000">0.093</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>286 (46)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.09 (0.93, 1.29)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.28</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>20 (54)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.42 (0.75, 2.68)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.28</font></td>
</tr> 
<tr><td colspan='18' style='font-weight: 900;'>Blood type</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;AB</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>317 (40)</td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>21 (32)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>21 (32)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>261 (42)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>8 (22)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'></td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;O</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>102 (13)</td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>10 (15)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.15 (0.54, 2.45)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#320000">0.72</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>10 (15)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.48 (0.70, 3.12)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.31</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>77 (12)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>0.70 (0.54, 0.91)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#FF0000">0.008</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>6 (16)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>2.35 (0.82, 6.74)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.11</font></td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>40 (5)</td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>3 (5)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.18 (0.35, 3.96)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#320000">0.79</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>3 (5)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>1.14 (0.34, 3.87)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#320000">0.83</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>32 (5)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>0.91 (0.59, 1.39)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#320000">0.65</font></td>
<td style='' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>3 (8)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'>3.04 (0.81, 11.43)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; text-align: center;'><font color="#650000">0.10</font></td>
</tr>
<tr>
<td style='border-bottom: 2px solid grey; text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>338 (42)</td>
<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>32 (48)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>0.99 (0.57, 1.73)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'><font color="#320000">0.98</font></td>
<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>32 (48)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>1.48 (0.85, 2.56)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'><font color="#650000">0.17</font></td>
<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>248 (40)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>0.59 (0.49, 0.70)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'><font color="#FF0000">&lt; 0.001</font></td>
<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>20 (54)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'>2.37 (1.04, 5.41)</td>
<td style='white-space: nowrap; padding: 0px 5px 0px; border-bottom: 2px solid grey; text-align: center;'><font color="#CC0000">0.040</font></td>
</tr>
</tbody>
</table>

### cuminc formula method


```r
## can use same formula as crr2
ci1 <- cuminc2(Surv(futime, event(censored) == death) ~ abo,
               data = transplant)

## but event indicator is not required
ci2 <- cuminc2(Surv(futime, event(censored)) ~ sex,
               data = transplant)
```

### cuminc2 summary methods


```r
summary(ci1)
```

```
## $est
##                   0     500    1000    1500    2000
## A death     0.00000 0.06191 0.06191      NA      NA
## B death     0.00000 0.08799 0.08799 0.12717 0.12717
## AB death    0.00000 0.07317      NA      NA      NA
## O death     0.00289 0.09106 0.09598 0.09598      NA
## A ltx       0.00000 0.81022 0.84174      NA      NA
## B ltx       0.00000 0.71905 0.77478 0.77478 0.77478
## AB ltx      0.02439 0.80488      NA      NA      NA
## O ltx       0.00000 0.69130 0.76563 0.78129      NA
## A withdraw  0.00000 0.02481 0.02481      NA      NA
## B withdraw  0.00000 0.05886 0.05886 0.05886 0.05886
## AB withdraw 0.00000 0.07317      NA      NA      NA
## O withdraw  0.00000 0.04413 0.06470 0.06470      NA
## 
## $var
##                   0     500    1000    1500    2000
## A death     0.00000 0.00018 0.00018      NA      NA
## B death     0.00000 0.00081 0.00081 0.00259 0.00259
## AB death    0.00000 0.00192      NA      NA      NA
## O death     0.00001 0.00024 0.00027 0.00027      NA
## A ltx       0.00000 0.00048 0.00046      NA      NA
## B ltx       0.00000 0.00202 0.00233 0.00233 0.00233
## AB ltx      0.00059 0.00426      NA      NA      NA
## O ltx       0.00000 0.00064 0.00058 0.00062      NA
## A withdraw  0.00000 0.00008 0.00008      NA      NA
## B withdraw  0.00000 0.00056 0.00056 0.00056 0.00056
## AB withdraw 0.00000 0.00177      NA      NA      NA
## O withdraw  0.00000 0.00012 0.00020 0.00020      NA
## 
## $events
##             0 500 1000 1500 2000
## A death     0  20   20   21   21
## B death     0   9    9   10   10
## AB death    0   3    3    3    3
## O death     1  31   32   32   32
## A ltx       0 262  269  269  269
## B ltx       0  74   77   77   77
## AB ltx      1  33   33   33   33
## O ltx       0 234  254  256  256
## A withdraw  0   8    8    8    8
## B withdraw  0   6    6    6    6
## AB withdraw 0   3    3    3    3
## O withdraw  0  15   20   20   20
## 
## $total_events
##     A death     B death    AB death     O death       A ltx       B ltx      AB ltx       O ltx 
##          21          10           3          32         269          78          33         256 
##  A withdraw  B withdraw AB withdraw  O withdraw 
##           8           6           3          20 
## 
## $total_groups
##   A   B  AB   O 
## 325 103  41 346 
## 
## $total_atrisk
##    0  500 1000 1500 2000 
##  811   94   20    2    1
```

```r
summary(ci1, times = 0:10 * 100)$est
```

```
##                   0     100     200     300     400     500     600     700     800     900    1000
## A death     0.00000 0.03390 0.04317 0.05562 0.06191 0.06191 0.06191 0.06191 0.06191 0.06191 0.06191
## B death     0.00000 0.04854 0.06796 0.07767 0.08799 0.08799 0.08799 0.08799 0.08799 0.08799 0.08799
## AB death    0.00000 0.04878 0.04878 0.04878 0.04878 0.07317      NA      NA      NA      NA      NA
## O death     0.00289 0.05250 0.07597 0.08185 0.08779 0.09106 0.09106 0.09106 0.09106 0.09106 0.09598
## A ltx       0.00000 0.48743 0.72831 0.77824 0.80346 0.81022 0.83623 0.84174 0.84174 0.84174 0.84174
## B ltx       0.00000 0.33010 0.50485 0.69903 0.71905 0.71905 0.74866 0.74866 0.74866 0.77478 0.77478
## AB ltx      0.02439 0.56098 0.60976 0.73171 0.80488 0.80488      NA      NA      NA      NA      NA
## O ltx       0.00000 0.23354 0.43015 0.54777 0.64341 0.69130 0.73601 0.75159 0.75580 0.76071 0.76563
## A withdraw  0.00000 0.00926 0.01853 0.02167 0.02481 0.02481 0.02481 0.02481 0.02481 0.02481 0.02481
## B withdraw  0.00000 0.02913 0.03883 0.04854 0.05886 0.05886 0.05886 0.05886 0.05886 0.05886 0.05886
## AB withdraw 0.00000 0.02439 0.07317 0.07317 0.07317 0.07317      NA      NA      NA      NA      NA
## O withdraw  0.00000 0.01167 0.02927 0.03810 0.04413 0.04413 0.04413 0.06470 0.06470 0.06470 0.06470
```

### cuminc plotting methods


```r
par(mfrow = c(2, 2))

# ciplot(ci2)
plot(ci2, add = TRUE) ## equivalently

plot(ci2, split = 'event', add = TRUE, wh.events = 'est')
```

![](cmprk2_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
## convenience wrapper
par(mfrow = c(2, 2))
ciplot_by(
  rhs = 'sex', time = 'futime', event = 'event',
  data = transplant, by = 'abo', xlim = c(0, 1500),
  events = FALSE, single = FALSE, events.total = 2100
)
```

![](cmprk2_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

### extras


```r
## pairwise gray tests
cuminc_pairs(ci1)$p.value
```

```
## $death
##        A     B    AB  O
## A     NA 1.000 1.000  1
## B  0.391    NA 1.000  1
## AB 0.788 0.776    NA  1
## O  0.201 0.877 0.712 NA
## 
## $ltx
##        A     B    AB     O
## A     NA 0.035 0.805 0.000
## B  0.007    NA 0.442 0.442
## AB 0.805 0.162    NA 0.051
## O  0.000 0.147 0.013    NA
## 
## $withdraw
##        A     B    AB     O
## A     NA 0.446 0.446 0.253
## B  0.097    NA 1.000 1.000
## AB 0.089 0.753    NA 1.000
## O  0.042 0.885 0.402    NA
```

```r
timepoints2(
  ci2, html = TRUE,
  htmlArgs = list(
    caption = 'cuminc estimates at specific time points (<code>cuminc::timepoints</code>).'
  )
)
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr><td colspan='6' style='text-align: left;'>
cuminc estimates at specific time points (<code>cuminc::timepoints</code>).</td></tr>
<tr>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>0</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>500</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1000</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1500</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2000</th>
</tr>
</thead>
<tbody>
<tr>
<td style='text-align: left;'>m death</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.002</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.091</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.097</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.104</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.104</td>
</tr>
<tr>
<td style='text-align: left;'>f death</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.000</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.063</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.063</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>-</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>m ltx</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.000</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.738</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.786</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.800</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.800</td>
</tr>
<tr>
<td style='text-align: left;'>f ltx</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.003</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.761</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.822</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>-</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>m withdraw</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.000</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.030</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.045</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.045</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; text-align: center;'>0.045</td>
</tr>
<tr>
<td style='border-bottom: 2px solid grey; text-align: left;'>f withdraw</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; border-bottom: 2px solid grey; text-align: center;'>0.000</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; border-bottom: 2px solid grey; text-align: center;'>0.052</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; border-bottom: 2px solid grey; text-align: center;'>0.056</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; border-bottom: 2px solid grey; text-align: center;'>-</td>
<td style='padding: 0px 5px 0px; white-space: nowrap; border-bottom: 2px solid grey; text-align: center;'>-</td>
</tr>
</tbody>
</table>

```r
AIC(cr1$`CRR: death`)
```

```
##      AIC 
## 867.6161
```

```r
sapply(cr1, BIC)
```

```
##    CRR: death.BIC      CRR: ltx.BIC CRR: withdraw.BIC 
##          878.5644         7494.9446          496.6589
```

```r
logLik(cr1$`CRR: death`)
```

```
## 'log Lik.' -428.8081 (df=5)
```

```r
deviance(cr1$`CRR: death`)
```

```
## Call:
## crr(transplant[, "futime"], transplant[, "event"], cov1 = model.matrix(~age + 
##     sex + abo, transplant)[, -1L, drop = FALSE], cencode = "censored", 
##     failcode = "death", variance = TRUE, cengroup = rep(1L, nrow(transplant)), 
##     gtol = 1e-06, maxiter = 10, init = c(0L, 0L, 0L, 0L, 0L))
## 
## Deviance = 7.15 on 5 df, 0.20984
```

```r
crrFits(cr1$`CRR: death`)
```

```
## Model selection table
## 
## 0: Null Model 
## 
## 1: Model 1 call:
## crr(transplant[, "futime"], transplant[, "event"], cov1 = model.matrix(~age + 
##     sex + abo, transplant)[, -1L, drop = FALSE], cencode = "censored", 
##     failcode = "death", variance = TRUE, cengroup = rep(1L, nrow(transplant)), 
##     gtol = 1e-06, maxiter = 10, init = c(0L, 0L, 0L, 0L, 0L))
## 
## 
##     n  loglik df k -2logLik -2logLik diff    AIC AIC diff    BIC BIC diff
## 0 797 -432.38  0 1   864.76        7.1484 864.76   0.0000 864.76      0.0
## 1 797 -428.81  5 6   857.62        0.0000 867.62   2.8516 878.56     13.8
```

```r
crrwald.test(cr1$`CRR: death`)
```

```
##               chi2 df          P
## Overall 8.11151956  5 0.15019572
## age     2.86948310  1 0.09027386
## sexf    2.82220518  1 0.09296860
## aboB    1.04828004  1 0.30590354
## aboAB   0.04546327  1 0.83115446
## aboO    1.91025511  1 0.16693492
```
