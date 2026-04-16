# Append matched datasets together

These functions are [`rbind()`](https://rdrr.io/r/base/cbind.html)
methods for objects resulting from calls to
[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
and
[`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md).
They function nearly identically to
[`rbind.data.frame()`](https://rdrr.io/r/base/cbind.html); see Details
for how they differ.

## Usage

``` r
# S3 method for class 'matchdata'
rbind(..., deparse.level = 1)

# S3 method for class 'getmatches'
rbind(..., deparse.level = 1)
```

## Arguments

- ...:

  Two or more `matchdata` or `getmatches` objects the output of calls to
  [`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md)
  and
  [`get_matches()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md),
  respectively. Supplied objects must either be all `matchdata` objects
  or all `getmatches` objects.

- deparse.level:

  Passed to [`rbind()`](https://rdrr.io/r/base/cbind.html).

## Value

An object of the same class as those supplied to it (i.e., a `matchdata`
object if `matchdata` objects are supplied and a `getmatches` object if
`getmatches` objects are supplied).
[`rbind()`](https://rdrr.io/r/base/cbind.html) is called on the objects
after adjusting the variables so that the appropriate method will be
dispatched corresponding to the class of the original data object.

## Details

[`rbind()`](https://rdrr.io/r/base/cbind.html) appends two or more
datasets row-wise. This can be useful when matching was performed
separately on subsets of the original data and they are to be combined
into a single dataset for effect estimation. Using the regular
`data.frame` method for [`rbind()`](https://rdrr.io/r/base/cbind.html)
would pose a problem, however; the `subclass` variable would have
repeated names across different datasets, even though units only belong
to the subclasses in their respective datasets. `rbind.matchdata()`
renames the subclasses so that the correct subclass membership is
maintained.

The supplied matched datasets must be generated from the same original
dataset, that is, having the same variables in it. The added components
(e.g., weights, subclass) can be named differently in different datasets
but will be changed to have the same name in the output.

`rbind.getmatches()` and `rbind.matchdata()` are identical.

## See also

[`match_data()`](https://kosukeimai.github.io/MatchIt/reference/match_data.md),
[`rbind()`](https://rdrr.io/r/base/cbind.html)

See `vignettes("estimating-effects")` for details on using
[`rbind()`](https://rdrr.io/r/base/cbind.html) for effect estimation
after subsetting the data.

## Author

Noah Greifer

## Examples

``` r
data("lalonde")

# Matching based on race subsets
m.out_b <- matchit(treat ~ age + educ + married +
                    nodegree + re74 + re75,
                  data = subset(lalonde, race == "black"))
#> Warning: Fewer control units than treated units; not all treated units will get
#> a match.
md_b <- match_data(m.out_b)

m.out_h <- matchit(treat ~ age + educ + married +
                    nodegree + re74 + re75,
                  data = subset(lalonde, race == "hispan"))
md_h <- match_data(m.out_h)

m.out_w <- matchit(treat ~ age + educ + married +
                    nodegree + re74 + re75,
                  data = subset(lalonde, race == "white"))
md_w <- match_data(m.out_w)

#Bind the datasets together
md_all <- rbind(md_b, md_h, md_w)

#Subclass conflicts are avoided
levels(md_all$subclass)
#>   [1] "1_1"  "1_2"  "1_3"  "1_4"  "1_5"  "1_6"  "1_7"  "1_8"  "1_9"  "1_10"
#>  [11] "1_11" "1_12" "1_13" "1_14" "1_15" "1_16" "1_17" "1_18" "1_19" "1_20"
#>  [21] "1_21" "1_22" "1_23" "1_24" "1_25" "1_26" "1_27" "1_28" "1_29" "1_30"
#>  [31] "1_31" "1_32" "1_33" "1_34" "1_35" "1_36" "1_37" "1_38" "1_39" "1_40"
#>  [41] "1_41" "1_42" "1_43" "1_44" "1_45" "1_46" "1_47" "1_48" "1_49" "1_50"
#>  [51] "1_51" "1_52" "1_53" "1_54" "1_55" "1_56" "1_57" "1_58" "1_59" "1_60"
#>  [61] "1_61" "1_62" "1_63" "1_64" "1_65" "1_66" "1_67" "1_68" "1_69" "1_70"
#>  [71] "1_71" "1_72" "1_73" "1_74" "1_75" "1_76" "1_77" "1_78" "1_79" "1_80"
#>  [81] "1_81" "1_82" "1_83" "1_84" "1_85" "1_86" "1_87" "2_1"  "2_2"  "2_3" 
#>  [91] "2_4"  "2_5"  "2_6"  "2_7"  "2_8"  "2_9"  "2_10" "2_11" "3_1"  "3_2" 
#> [101] "3_3"  "3_4"  "3_5"  "3_6"  "3_7"  "3_8"  "3_9"  "3_10" "3_11" "3_12"
#> [111] "3_13" "3_14" "3_15" "3_16" "3_17" "3_18"
```
