![enter image description here][1]

#Implementing tests#


## Chi square equidistribution test ##
The first test we will examine is a relatively simple, ubiquitous test: the chi square test for independence. In this approach we will use the chi square test to measure how the observed values differ from expected values. 
###Zeros and ones###
The expected value distribution of zeros and ones is very simple. There is a 1/2 probability that a variate will be 0 and a 1/2 probability that a variate will be 1. Therefore, we can set up a statistic like the following:

    chiSquareTestStatIntegers[seq_] := (Total[seq] - 1/2 Length[seq])^2/(
       1/2 Length[
         seq])(*observed number of ones*)+ ((Length[seq] - Total[seq]) - 
         1/2 Length[seq])^2/(1/2 Length[seq])(*observed number of zeros*);

We know that this statistic will follow a chi square distribution with d.f. = 1. Let's confirm:

    Show[Histogram[
      Table[chiSquareTestStatIntegers[RandomInteger[{0, 1}, 10000]], 
       100000], Automatic, "PDF"], 
     Plot[PDF[ChiSquareDistribution[1], x], {x, 0, 10}, 
      PlotRange -> {0, 2}], Sequence[
     PlotLabel -> "Sampling distribtion of test statistic vs. chi square \
    distribution", Frame -> True, FrameLabel -> {"Test statistic", "PDF"}]
     ]

![enter image description here][2]

Therefore, we can construct a p-value as follows:

    chiSquareTestPValueIntegers[seq_] := 
     1. - CDF[ChiSquareDistribution[1], chiSquareTestStatIntegers[seq]]

We observe that the rate of type-I error behaves as expected for a well-functioning test. This, of course, assumes that the pseudorandom number generator that is `RandomInteger` is truly "random"—whatever that may mean.

### Random reals between zero and one ###

Again, with the chi square test for independence, we must compare observed values to expected values. For continuous distributions, we must choose bins in which we have expected values. For our purposes, we will choose bins of size 1/n where n is the length of the sequence. Using this bin specification, we know that the expected number of numbers we should observe in each bin should be 1. With this information, we can create a test statistic:

    chiSquareTestStatReals[seq_] := 
     Total[(BinCounts[seq, 1/Length[seq]] - 1)^2/
      1](*note we use the listable atribute of squaring and basic \
    arithmetic*)

This distribution should follow the chi square distribution with d.f. = n -2 where n is the length of the sequence. Let's confirm:

    Show[Histogram[
      Table[chiSquareTestStatReals[RandomReal[{0, 1}, 1000]], 10000], 
      Automatic, "PDF"], 
     Plot[PDF[ChiSquareDistribution[998], x], {x, 0, 1500}, 
      PlotRange -> {1.*^-6, All}], Sequence[
     PlotLabel -> "Sampling distribtion of test statistic vs. chi square \
    distribution", Frame -> True, FrameLabel -> {"Test statistic", "PDF"}]
     ]

![enter image description here][3]
## The Kolmogorov–Smirnov distribution fit test##

The Kolmogorov–Smirnov test is a powerful distribution fit test that compares an empirical CDF to a known CDF. The Kolmogorov–Smirnov test is relatively simple to formulate, but fortunately, it is already implemented in the Wolfram Language. We can then apply it to test for equidistribution of a sequence:

    sequence = RandomReal[{0, 1}, 10000]; 
    KolmogorovSmirnovTest[sequence, UniformDistribution[]]


This expression returns a p–value.

## Serial test ##
The serial test examines the frequency of the following overlapping pairs within the overall sequence: 

    {0, 0}      {0, 1}      {1, 0}      {1, 1}

For example, the sequence 

    {1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0,
       1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 
      1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0,
       0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 
      0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0}

has the associated counts

    {16, 23, 24, 21}

where the order is that of the pairs list above. We can create a chi square–like statistic that measures the deviation from expected frequency:

![enter image description here][4]

where n is the sequence length and µ_freq is the expected frequency with which the the pair occurs. These frequencies are computed empirically: 

    mcts[n_, samplesize_] := 
      Mean[Table[
        N@SequenceCount[RandomInteger[{0, 1}, n], #] & /@ 
         Tuples[{0, 1}, 2], samplesize]];

We can observe how these expected frequencies vary with sequence length:

    {pts00, pts01, pts10, pts11} = 
      ParallelTable[{{n, n, n, n}, mcts[n, 10000]/n}\[Transpose], {n, 500,
          10000, 500}]\[Transpose];
    ListPlot[{pts00, pts01, pts10, pts11}, Frame -> True, 
     FrameLabel -> {"Sequence length", "Mean pair counts ratio"}]

![ ][5]

This figure displays how the expected frequencies are invariant under sequence length variation. After discovering this, we can come up with one set of constants that will work for any sequence length:

    In[]:= $mctsrat = mcts[1000, 100000]/1000
    Out[]:={0.166564, 0.249761, 0.249774, 0.166603}

And with this in mind we can finally construct a test statistic:

    SerialTestStat[seq_] := 
      Total[((SequenceCount[seq, #] & /@ Tuples[{0, 1}, 2]) - 
         Length[seq] $mctsrat)^2/(Length[seq] $mctsrat)];

As we chose a chi square\[Dash]like statistic, we can anticipate that the test statistic will follow a gamma distribution.

    serialsamplingdist[n_, samplesize_] := 
     ParallelTable[SerialTestStat[RandomInteger[{0, 1}, n]], samplesize]
    serialdistsample1000 = serialsamplingdist[1000, 100000];
    Histogram[serialdistsample1000, {.1}, "PDF"]

![enter image description here][6]

    In[]:= serialfitdist1000 = 
     FindDistribution[serialdistsample1000, 
      TargetFunctions -> {GammaDistribution}]
    Out[]:= GammaDistribution[1.12212, 1.53172] 
.

    Show[Histogram[serialdistsample1000, {.1}, "PDF"], 
     Plot[PDF[serialfitdist1000, x], {x, 0, 10}]]

![enter image description here][7]

We can see that this distribution fits relatively well. Let's now investigate how the distribution parameters \[Alpha] and \[Beta] vary with sequence length.

    \[Alpha]\[Beta]params[n_, samplesize_] := 
     FindDistributionParameters[serialsamplingdist[n, samplesize], 
      GammaDistribution[\[Alpha], \[Beta]]]
    paramvslength = 
     Table[{{n, 
          n}, {\[Alpha], \[Beta]}}\[Transpose] /. \[Alpha]\[Beta]params[n,
         10000], {n, 100, 10000, 1000}]
    Row@{ListPlot[paramvslength\[Transpose][[1]], Sequence[
       ImageSize -> Medium, PlotLabel -> "\[Alpha]", 
        PlotRange -> {0, All}]], 
      ListPlot[paramvslength\[Transpose][[2]], Sequence[
       ImageSize -> Medium, PlotLabel -> "\[Beta]", 
        PlotRange -> {0, All}]]}

![ ][8]

We can see that the parameters hardly vary and do not vary consistently with list length. We can then find parameters we like.

    In[]:= $serialparams = \[Alpha]\[Beta]params[10000, 10000]
    Out[]= {\[Alpha] -> 1.10793, \[Beta] -> 1.53947}

    In[]:= $serialfitdist = 
     GammaDistribution[\[Alpha], \[Beta]] /. $serialparams
    Out[]= GammaDistribution[1.10793, 1.53947]

Our p–value is then of the form:

    1 - CDF[$serialfitdist, SerialTestStat[sequence]]

###Discussion###
This method is powerful as it examines the frequency and ordering of pairs. However, this method is contingent on one large assumption\[LongDash]that the Wolfram Language pseudorandom number generator is truly pseudorandom. Or rather, that it is pseudorandom enough. A Monte Carlo approach is often very powerful for these sorts of statistics; however, it is important to be weary of such assumptions.

## Runs tests ##

In these tests, we will investigate runs of both 0s and 1s and random reals between 0 and 1. A run is a sequence that fits a particular pattern until it is reset. For random reals, we will use increasing runs. For sequences of 0s and 1s, we obtain the runs by simply splitting the list:

###How long are runs in a sequence?###

The first statistic we will construct will measure the length of runs. Our approach will again utilize a Monte Carlo method. To measure runs:

    runsup[seq_] := Split[seq, #1 < #2 &];
In order to count lengths of runs:

    runLengthDist[seq_] := Length /@ runsup[seq]
The distribution of run lengths within a random sequence:

    ListPlot[Counts@runLengthDist@RandomReal[{0, 1}, 1000000], 
     Filling -> Axis]

![enter image description here][9]
 
We can see that the run lengths quickly converge to a constant distribution. With this we can then calculate a constant–the mean run length:

    In[]:= $meanRunLength = 
     Mean[ParallelMap[ N@Mean[runLengthDist[#]] &, 
       RandomReal[{0, 1}, {10000, 5000}]]] 
    Out[]:= 2.00004 

We will then construct another chi\[Dash]square like test statistic:

    RunLengthStat[seq_] := 
      Total[(runLengthDist[seq] - $meanRunLength)^2/$meanRunLength];

This statistic turns out to be normally distributed with a mean and standard deviation associated with the sequence length. Putting all of this together, we can calculate the two–tailed p–value:

    2 Min[{1 - CDF[##], CDF[##]}] &[
     NormalDistribution[runlength\[Mu][Length[sequence]], 
      runlength\[Sigma][Length[sequence]]], RunLengthStat[sequence]];

###How many runs do we expect a sequence to have?###
The next statistic we create will count the total amount of runs a sequence has. The statistic is then: 

    RunCountStat = Length[Split[#, (#1 < #2 &)]] &;

This distribution follows a normal distribution:

    Manipulate[
     ProbabilityPlot[
      RunCountStat /@ RandomReal[{0, 1}, {1000, seqleng}]], {seqleng, 
      5000}]

![enter image description here][10]

We can calculate and model the distribution using a Monte–Carlo method to create a two–tailed p–value:

    2 Min[{1 - CDF[##], CDF[##]}] &[
     NormalDistribution[runcount\[Mu][Length[sequence]], 
      runcount\[Sigma][Length[sequence]]], RunCountStat[sequence]] 
where

    runcount\[Mu][seqleng_] = seqleng/2
    runcount\[Sigma][seqleng_] = 0.303833 seqleng^0.493785 

###How many runs does a sequence of zeros and ones have?###

This test is taken entirely from [an NIST publication][11]. It only assesses sequences of 0s and 1s. Additionally, it is theoretically derived. The test statistic the publication gives is

    NISTRunsstat = Length@Split@# &;

And the publication assigns the associated p\[Dash]value.

    NISTRunspval = 
      Erfc[Abs[Length@Split@#1 - 2 Length[#1] #2 (1 - #2)]/(
          2. Sqrt[2 Length[#1]] #2 (1 - #2))] &[#, Total[#]/Length[#]] &;

##Spectral tests##

The spectral test will use the discrete cosine transform to assess the randomness of data. The transformed data should be homogeneous over frequency space and, overall, should be normally distributed. The transformed data should be normally distributed because the transformed data is obtained by:

![enter image description here][12]

where u[r] is the rth random number in the sequence, s is the frequency parameter, and n is the length of the sequence. This can be re–expressed in the form:

![enter image description here][13]

where U is the sequence, a[s] is a list of constants, and \[CenterDot] is the dot product. This form tells us that each individual point is a sum of random variables. If these variables are truly independently and identically distributed, they should follow the normal distribution:

    sequence = RandomReal[{0, 1}, 10000];
    ListPlot[Rest@FourierDCT@sequence, PlotRange -> All]
    Histogram[Rest@FourierDCT@sequence, Automatic, "PDF"] 

![enter image description here][14]
![enter image description here][15]

Our task then is to simply perform a Kolmogorov\[Dash]Smirnov test for normality on the transformed data:

    KolmogorovSmirnovTest[Rest@FourierDCT@sequence, Automatic, "PValue"]

##Sequential sums tests##

This test is entirely developed by the NIST publication. The test statistic for this test is based on the idea of generating a random walk from a list of zeros and ones. First, we transform the sequence and then accumulate the values.

    sequence = RandomInteger[{0, 1}, 1000];
    ListLinePlot[Accumulate[2 sequence - 1]]

![enter image description here][16]

For the statistic, we will look at the maximum distance away from zero that the walk obatins:

    CUSUMmaxstat[seq_] := Max@Abs@Accumulate[2 seq - 1];

The stationary distribution is approximately normal. The publication uses this property to derive the following formula for the p–value:

    \[CapitalPhi][x_] = CDF[NormalDistribution[], x];
    CUSUMmaxpval = 
      Function[{z, n}, 
         1. - Sum[\[CapitalPhi][(z (4 k + 1))/Sqrt[n]] - \[CapitalPhi][(
             z (4 k - 1))/Sqrt[n]], {k, Floor[(-n/z + 1)/4], (n/z - 1)/4, 
            1}] + Sum[\[CapitalPhi][(z (4 k + 3))/Sqrt[
             n]] - \[CapitalPhi][(z (4 k + 1))/Sqrt[n]], {k, 
            Floor[(-n/z - 3)/4], (n/z - 1)/4, 1}]][CUSUMmaxstat[#], 
        Length[#]] &;

##Arcsine test##
This test utilizes the arcsine law as described in Lorek, Los, Gotfryd, and Zagórski's publication. The arcsine law states that given a variable X—the percentage of a random walk generated by a sequence of random reals spends above (or equivalently below) 0—X will be distributed according to a particular 'arcsine distribution':

    Table[Length@
       Cases[Accumulate[2 # - 1] &@RandomReal[{0, 1}, 1000], x_ /; x > 0]/
      1000, 10000];
    Show[Histogram[%, 100, "CDF"], 
     Plot[2/\[Pi] ArcSin[Sqrt[x]], {x, 0, 1}]]

![enter image description here][17]

The formula for the CDF of the distribution is 

![enter image description here][18]

Our statistic is then the percentage of time a random walk generated from our sequence spends above 0:

    Arcsinelawstat[seq_] := 
      Length@Cases[Accumulate[2 # - 1] &@seq, x_ /; x > 0]/Length[seq];

And we create a two\[Dash]tailed p\[Dash]value:

    Arcsinelawpval[stat_] := 
      2 Min@{1 - #, #} &[2/\[Pi] ArcSin[Sqrt[stat]]];

  [1]: https://community.wolfram.com//c/portal/getImageAttachment?filename=4860ScreenShot2019-07-07at3.54.16PM.png&userId=1705502
  [2]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled.png&userId=1705502
  [3]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled1.png&userId=1705502
  [4]: https://community.wolfram.com//c/portal/getImageAttachment?filename=ScreenShot2019-07-07at4.14.24PM.png&userId=1705502
  [5]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled2.png&userId=1705502
  [6]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled4.png&userId=1705502
  [7]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled5.png&userId=1705502
  [8]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled6.png&userId=1705502
  [9]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled7.png&userId=1705502
  [10]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled11.png&userId=1705502
  [11]: https://csrc.nist.gov/publications/detail/sp/800-22/rev-1a/final
  [12]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled12.png&userId=1705502
  [13]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled13.png&userId=1705502
  [14]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled14.png&userId=1705502
  [15]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled15.png&userId=1705502
  [16]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled16.png&userId=1705502
  [17]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled17.png&userId=1705502
  [18]: https://community.wolfram.com//c/portal/getImageAttachment?filename=Untitled18.png&userId=1705502